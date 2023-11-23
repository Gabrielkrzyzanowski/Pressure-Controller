#!/usr/bin/env python3  

import itertools
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import spe_loader as sl 
import numpy as np 
import random
from lmfit import models 
import time
from epics import caput,caget 
import os 
from threading import Thread
import argparse
from scipy.signal import find_peaks 

file_path='/ibira/lnls/beamlines/ema/apps/LightField/Experiments_Files'

prefix = 'EMA:B:PACE01:'
pv_setpoint = prefix + 'Setpoint'
pv_pressure_rbv = prefix + 'Pressure_RBV' 
pv_acquire_mode='EMA:13LF1:cam1:Acquire'
pv_file_number='EMA:13LF1:cam1:FileNumber'

time1=time.time()

def Model(spec): # Chris Ostrouchov's function to create a robust Model for XRD data fitting
    composite_model=None 
    params = None
    x=spec['x']
    y=spec['y']
    
    x_max=np.max(x)
    x_min=np.min(x)
    x_range=x_max-x_min
    y_max=np.max(y) 

    for i,function in enumerate(spec['model']):
        prefix=f'model{i}_'
        model=getattr(models,function['type'])(prefix=prefix)
            
        if function['type'] in ['GaussianModel', 'LorentzianModel', 'VoigtModel','LinearModel']:
            
            model.set_param_hint('sigma', min=1e-6, max=x_range)
            model.set_param_hint('center', min=x_min, max=x_max)
            model.set_param_hint('height', min=1e-6, max=1.1*y_max)
            model.set_param_hint('amplitude', min=1e-6)

            default_params = {
                prefix+'center': x_min + x_range * random.random(),
                prefix+'height': y_max * random.random(),
                prefix+'sigma': x_range * random.random()
            } 
        else:
            raise NotImplemented(f'model {function["type"]} not implemented yet')
        if 'help' in function:  # allow override of settings in parameter
            for param, options in function['help'].items():
                model.set_param_hint(param, **options)

        model_params = model.make_params(**default_params, **function.get('params', {}))

        if params is None:
            params = model_params
        else:
            params.update(model_params)
        if composite_model is None:
            composite_model = model
        else:
            composite_model = composite_model + model

    return composite_model,params

def PeakFinder(wave_lenght,ruby_data): 
    """
    Find and fit the optimal peaks in the given data.

    Parameters:
    - wave_length (array): The wavelength data.
    - ruby_data (array): The corresponding data values.

    Returns:
    - result (lmfit.ModelResult): The result of the peak fitting.
    - peaks (array): The indices of the detected peaks.
    """

    peak_parameters={'prominence':1, 'width':3, 'height':1500}
    peaks, properties = find_peaks(ruby_data,
                                   height=peak_parameters['height'],
                                   prominence=peak_parameters['prominence'],
                                   width=peak_parameters['width'])  
    
    while len(peaks)>2: 
        peak_parameters['width'] += 0.5
        peaks, properties = find_peaks(ruby_data,
                                       height=peak_parameters['height'],
                                       prominence=peak_parameters['prominence'],
                                       width=peak_parameters['width']) 
    while len(peaks)<2: 
        peak_parameters['width'] -= 0.5
        peaks, properties = find_peaks(ruby_data,
                                       height=peak_parameters['height'],
                                       prominence=peak_parameters['prominence'],
                                       width=peak_parameters['width']) 
        
    spec = {
    'x': wave_lenght,
    'y': ruby_data,
    'model': [
        {
        'type': 'VoigtModel' ,
        'params':{'center': wave_lenght[peaks[0]], 'height': properties['peak_heights'][0], 'sigma': 0.25},
        'help': {'center': {'min': wave_lenght[peaks[0]] - peak_parameters['width']/2, 
                            'max': wave_lenght[peaks[0]] + peak_parameters['width']/2}}
        },
        {
        'type': 'VoigtModel' ,
        'params':{'center': wave_lenght[peaks[1]], 'height': properties['peak_heights'][1], 'sigma': 0.25},
        'help': {'center': {'min': wave_lenght[peaks[1]] - peak_parameters['width']/2, 
                            'max': wave_lenght[peaks[1]] + peak_parameters['width']/2}}
        }, 
        {'type': 'LinearModel'}
    ]
    } 

    model, params = Model(spec)
    result = fit_model(ruby_data,  wave_lenght, model, params)

    return result, peaks 

def fit_model(ruby_data,  wave_lenght, model, params, methods=['least_squares', 'powell'], max_nfev_values=[1000, 2000], r_square_threshold=0.95):
    for method in methods:
        for max_nfev in max_nfev_values:
            try:
                result = model.fit(ruby_data, x= wave_lenght, **params, max_nfev=max_nfev, method=method)
                r_square = 1 - result.residual.var() / np.var(ruby_data)
                if r_square > r_square_threshold:
                    return result
            except Exception as e:
                print(f"An error occurred: {e}")
    
    return None

def CalculatePressure(WaveLenght,RubySpec): 
    
    A = 1870
    B = 5.63
    lmbd0 = 694.26 

    try: 
        result, peaks = PeakFinder(WaveLenght,RubySpec)
        lmbd = round(WaveLenght[result.best_fit.argmax()],4)

        return A*((lmbd - lmbd0)/lmbd0)*(1 + B*((lmbd - lmbd0)/lmbd0))
    except: 
        return None 

class DACPressureValueThread(Thread):
    
    # constructor
    def __init__(self):
        # execute the base constructor
        Thread.__init__(self)
        # set a default value
        self.PressureValue = None
 
    # function executed in a new thread
    def run(self):
        global i
        ruby_intesity=[]
        caput(pv_acquire_mode,1,wait=True)
        
        experiment_files=os.listdir(file_path)
        experiment_files.sort()

        
        if len(experiment_files)<=3:

            while len(experiment_files)<3: 

                caput(pv_acquire_mode,1,wait=True) 

                experiment_files=os.listdir(file_path)
                experiment_files.sort()

        ruby_data=sl.SpeFile(f'{file_path}/{experiment_files[0]}')          

        ruby_intesity.append(ruby_data.data[0][0][0]) 
        wave_lenght = ruby_data.wavelength

        RubySpec = ruby_intesity[0]
        WaveLenght = np.resize(wave_lenght, len(RubySpec))
        
        WaveLenght_ROI=WaveLenght[np.where((WaveLenght > 690) & (WaveLenght < 720))] 
        RubySpec_ROI= RubySpec[np.where((WaveLenght > 690) & (WaveLenght < 720))] 

        os.remove(file_path +'/'+ experiment_files[0]) 

        if experiment_files[-1]=='test_970.spe': 
            caput(pv_file_number,1,wait=True) 
            caput(pv_acquire_mode,1,wait=True)

            for f in os.listdir(file_path):
                if f=='test_001_spe': pass
                else: os.remove(os.path.join(file_path, f))
        
        result= CalculatePressure(WaveLenght_ROI,RubySpec_ROI) #DAC Pressure Value   
        
        self.PressureValue = result if result is not None else self.PressureValue 

class PIDClass:
    
    def __init__(
            self,
            PressureSetPoint: float = None,
            Kp: float = None,
            Ki: float = None, 
            Kd: float= None,
            PressureAlarm: float=None,
            FilePath: str = None
            ) -> None: 
        
        super(PIDClass,self).__init__()
    
        self.PressureSetPoint=PressureSetPoint

        self.Kp=Kp if Kp!=None else 0.4239883 #Calibrated Value 
        self.Ki=Ki if Ki!=None else 0.0065433 #Calibrated Value 
        self.Kd=Kd if Kd!=None else 0.0000000 #Calibrated Value 

        self.PressureAlarm=PressureAlarm

        self.FilePath=FilePath 
        self.RealMbnPressure=caget(pv_pressure_rbv)
        self.mbn_threshold=caget(pv_pressure_rbv)
        self.step,self.controller_output ,self.pressure_ruby, self.mbn_pressure, self.experiment_time= [], [], [], [], []
        self.SetupAnimation()

    def graph_init(self):
        
        self.ax0.set_ylim(-1.0, 1.1)
        self.ax0.set_xlim(0, 1)

        self.ax1.set_ylim(-1.0, 1.1)
        self.ax1.set_xlim(0, 1) 

        self.ax1.set_ylim(-1.0, 1.1)
        self.ax1.set_xlim(0, 1)
        
        del self.pressure_ruby[:]
        del self.controller_output[:] 
        del self.mbn_pressure[:]
        del self.step[:] 

        self.line1.set_data(self.step,self.pressure_ruby)
        self.line2.set_data(self.step, self.controller_output)
        self.line3.set_data(self.step, self.mbn_pressure) 

        return self.line1,self.line2,self.line3,

    def SetupAnimation(self): 

        self.fig=plt.figure(figsize=(7,10))
        
        self.ax0,self.ax1,self.ax2 = self.fig.add_subplot(311),self.fig.add_subplot(312),self.fig.add_subplot(313)

        self.ax0.grid()
        self.ax1.grid()
        self.ax2.grid() 

        self.ax0.title.set_text('Ruby Pressure')
        self.ax0.set_ylabel('GPa')

        self.ax1.title.set_text('Controller Output')
        self.ax1.set_ylabel('bar')

        self.ax2.title.set_text('Membrane Pressure')
        self.ax2.set_ylabel('bar')

        self.ax2.set_xlabel('Step')
        
        self.line1,= self.ax0.plot([],[])
        self.line2,= self.ax1.plot([],[])
        self.line3,= self.ax2.plot([],[])

        self.fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.3) 

        self.fig.suptitle(f'PID Calibration Experiment \n Parameters: kp={self.Kp}, ki={self.Ki}, kd={self.Kd}, Pressure Set Point={self.PressureSetPoint} GPa')

        ani = animation.FuncAnimation(self.fig, self.run, self.generator, interval=1000, init_func=self.graph_init,
                        save_count=1000000)

        ani.save('Animation.png') #This line will save the animation after 10000000 iteration. In practice the script will run from a long time and you can stop it at any time

    def run(self,data): 

        global time1
        t, y1, y2, y3= data
        
        time2=time.time() 
        experiment_time=round(time2-time1,3)  

        self.step.append(t)
        self.pressure_ruby.append(y1)
        self.controller_output.append(y2) 
        self.mbn_pressure.append(y3)
        self.experiment_time.append(experiment_time)

        np.savetxt(f'Tests/{self.FilePath}.csv',  [p for p in zip(self.step, self.controller_output, self.pressure_ruby,self.mbn_pressure,self.experiment_time)],delimiter=',', fmt='%5.4f',
                header='{0:^5s},{1:^7s},{2:^9s},{3:^11s},{4:^13s}'.format('Step','Controler Setpoint (bar)','DAC Pressure (GPa)','Membrane Pressure (bar)','Experiment Time (s)',))
        
        print(f'-----------\nStep: {self.step[-1]}') 
        print(f'Controller Output: {self.controller_output[-1]}')
        print(f'Mbn Real Pressure: {self.mbn_pressure[-1]}')
        print(f'Tempo de experimento: {experiment_time}') 
        print(f'DAC Pressure: {self.pressure_ruby[-1]}\n-----------') 
        

        xmin1, xmax1 = self.ax0.get_xlim() 
        xmin2, xmax2 = self.ax1.get_xlim()
        xmin3, xmax3 = self.ax2.get_xlim() 

        ymin1, ymax1 = self.ax0.get_ylim() 
        ymin2, ymax2 = self.ax1.get_ylim()
        ymin3, ymax3 = self.ax2.get_ylim()

        #Adjusting x limits
        if t >= xmax1:
            self.ax0.set_xlim(xmin1, 2*xmax1)
            self.ax0.figure.canvas.draw()
        
        if t >= xmax2:
            self.ax1.set_xlim(xmin2, 2*xmax2)
            self.ax1.figure.canvas.draw()

        if t >= xmax3:
            self.ax2.set_xlim(xmin3, 2*xmax3)
            self.ax2.figure.canvas.draw() 

        #Adjusting y limits
        if y1 >= ymax1:
            self.ax0.set_ylim(ymin1, 1.5*ymax1)
            self.ax0.figure.canvas.draw()
        
        if y2 >= ymax2:
            self.ax1.set_ylim(ymin2, 1.5*ymax2)
            self.ax1.figure.canvas.draw()

        if y3 >= ymax3:
            self.ax2.set_ylim(ymin3, 1.5*ymax3)
            self.ax2.figure.canvas.draw()

        
        self.line1.set_data(self.step,self.pressure_ruby)
        self.line2.set_data(self.step, self.controller_output)
        self.line3.set_data(self.step, self.mbn_pressure)  

        return self.line1,self.line2,self.line3,         

    def generator(self):  
        caput('EMA:B:PACE01:Control', 1)
        i=0

        self.PressureRBVThread= DACPressureValueThread() 
        self.PressureRBVThread.start()
        self.PressureRBVThread.join()
        self.PressureRBV=round(self.PressureRBVThread.PressureValue,4) 

        controler = self.LoopPID(self.Kp,self.Ki,self.Kd) 
        controler.send(None)
        self.SetPressure = controler.send((i,self.PressureRBV, self.PressureSetPoint))
        
        while True:
             
            self.PressureRBVThread= DACPressureValueThread() 
            self.PressureRBVThread.start()
            self.PressureRBVThread.join()
            self.PressureRBV=round(self.PressureRBVThread.PressureValue,4)

            controler = self.LoopPID(self.Kp,self.Ki,self.Kd) 
            controler.send(None)

            self.SetPressure = controler.send((i,self.PressureRBV, self.PressureSetPoint))
            
            if self.PressureRBV >= self.PressureAlarm: 
                print('ALARM PRESSURE EXCEED!!!!!! \nClosing the script...')
                print('Pressure Value: ',self.Pressure_RBV)
                #raise SystemExit()
            else:   
                caput(pv_setpoint, self.SetPressure, wait=True)
 
            self.RealMbnPressure= caget(pv_pressure_rbv)  
            
            i+=1  

            yield i, self.PressureRBV, self.SetPressure, self.RealMbnPressure

    def LoopPID(self,Kp,Ki,Kd):
        #Initializing
        e_prev=0 
        t_prev=-100
        SetPressure=0
        P,I,D=0,0,0

        while True:
            # yield SetPressure, wait for new t, PressureRBV and PressureSetPoint
            t , PressureRBV , PressureSetPoint = yield SetPressure

            error = PressureSetPoint - PressureRBV #Error 

            P += Kp*error
            I += Ki*error*(t-t_prev)
            D += Kd*(error-e_prev)/(t-t_prev)

            SetPressure = self.mbn_threshold+P+I+D

            #update stored data for next iteration
            e_prev=error 
            t_prev=t 

if __name__ == '__main__': 
        
    parser = argparse.ArgumentParser(description = 'Script to test the PID Controler.')

    #Positional Arguments
    parser.add_argument('pressure_set_point', type=float, help='The pressure setpoint (in GPa) sent to the PID')
    parser.add_argument('Kp', type=float, help='Proportional gain')
    parser.add_argument('Ki', type=float, help='Integral gain')
    parser.add_argument('Kd', type=float, help='Derivative gain')
    parser.add_argument('pressure_alarm', type=int, help='The pressure (in GPa) that must not be exceed')
    parser.add_argument('file_name', type=str, help='The name of the file where the data will be saved ')

    args = parser.parse_args()

    PressureSetPoint=args.pressure_set_point
    Kp=args.Kp
    Ki=args.Ki
    Kd=args.Kd
    PressureAlarm=args.pressure_alarm
    FileName=args.file_name

    PIDClass(PressureSetPoint, Kp, Ki, Kd, PressureAlarm, FileName)
