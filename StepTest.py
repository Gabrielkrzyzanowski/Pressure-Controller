#!/usr/bin/env python3  

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import spe_loader as sl 

import numpy as np 
import time
from epics import caput,caget 
import os 
from threading import Thread 
import argparse 

from lmfit import models 
import random
from scipy.signal import find_peaks 

time1=time.time()

file_path='/ibira/lnls/beamlines/ema/apps/LightField/Experiments_Files'
ruby_intesity=[]

prefix = 'EMA:B:PACE01:'
pv_setpoint = prefix + 'Setpoint'
pv_setpoint_rbv = pv_setpoint + '_RBV'
pv_slew = prefix + 'Slew'
pv_pressure_rbv = prefix + 'Pressure_RBV'
pv_exposition_time=prefix +'something' 
pv_acquire_mode='EMA:13LF1:cam1:Acquire'
pv_file_number='EMA:13LF1:cam1:FileNumber'


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

def PeakFinder(wave_lenght,ruby_data): # Creating a composite model to fit our data 
    

    peak_parameters={'prominence':1, 'width':3, 'height':2000}
    peaks, properties = find_peaks(ruby_data,
                                   height=peak_parameters['height'],prominence=peak_parameters['prominence'],width=peak_parameters['width'])  
    
    while len(peaks)>2: 
        peak_parameters['width'] += 0.5
        peaks, properties = find_peaks(ruby_data,
                                height=peak_parameters['height'],prominence=peak_parameters['prominence'],width=peak_parameters['width']) 
    while len(peaks)<2: 
        peak_parameters['width'] -= 0.5
        peaks, properties = find_peaks(ruby_data,
                                height=peak_parameters['height'],prominence=peak_parameters['prominence'],width=peak_parameters['width']) 
        
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
    result = model.fit(ruby_data, x=wave_lenght, **params, max_nfev=500, method='least_squares') 
    
    try: 
        r_square = (1 - result.residual.var() / np.var(ruby_data)) 
        if r_square>0.95: 
            return result, peaks  
        else:
            result = model.fit(ruby_data, x=wave_lenght, **params, max_nfev=500, method='powell') 
            try: 
                r_square = (1 - result.residual.var() / np.var(ruby_data)) 
                return result, peaks if r_square>0.95 else None
            except: 
                result = model.fit(ruby_data, x=wave_lenght, **params, max_nfev=1000, method='powell')
                return result, peaks if (1 - result.residual.var() / np.var(ruby_data))>0.95 else None
    
    except:
        result = model.fit(ruby_data, x=wave_lenght, **params, max_nfev=500, method='powell') 
        try: 
            r_square = (1 - result.residual.var() / np.var(ruby_data)) 
            return result, peaks if r_square>0.95 else None
        except:
            result = model.fit(ruby_data, x=wave_lenght, **params, max_nfev=1000, method='powell')
            return result, peaks if (1 - result.residual.var() / np.var(ruby_data))>0.95 else None
    


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
                time.sleep(0.2) 
                experiment_files=os.listdir(file_path)
                experiment_files.sort()

        ruby_data=sl.SpeFile(f'{file_path}/{experiment_files[0]}')          

        ruby_intesity.append(ruby_data.data[0][0][0]) 
        wave_lenght = ruby_data.wavelength

        RubySpec = ruby_intesity[0]
        WaveLenght = np.resize(wave_lenght, len(RubySpec))
        
        
        #np.savetxt(f'RubySpecs/DAC600_StepTest3_01/RubySpec{i}.csv',  [p for p in zip(RubySpec, WaveLenght)],delimiter=',', fmt='%5.4f',
        #    header='{0:^5s},{1:^7s}'.format('Ruby Intensity','Wave Length (nm)',))
        
        
        os.remove(file_path +'/'+ experiment_files[0]) 

        if experiment_files[-1]=='test_970.spe': 
            caput(pv_file_number,1,wait=True) 
            caput(pv_acquire_mode,1,wait=True)
            time.sleep(0.2) 
            for f in os.listdir(file_path):
                if f=='test_001_spe': pass
                else: os.remove(os.path.join(file_path, f))
        
        result= CalculatePressure(WaveLenght,RubySpec) #DAC Pressure Value  
        
        self.PressureValue = result if result != None else self.PressureValue 
           

class PIDClass:
    
    def __init__(
            self,
            PressureSetPoint: float = None,
            CurrentMbnPressure: float = None,
            FilePath: str = None
            ) -> None: 
        
        super(PIDClass,self).__init__()
    
        self.PressureSetPoint=PressureSetPoint
        self.CurrentMbnPressure=CurrentMbnPressure
        self.FilePath=FilePath 
        self.RealMbnPressure=self.CurrentMbnPressure  

        self.step, self.pressure_mbn ,self.pressure_ruby, self.mbn_real_pressure, self.experiment_time= [], [], [], [], []

        self.SetupAnimation()

    def graph_init(self):
        
        self.ax0.set_ylim(-1.0, 1.1)
        self.ax0.set_xlim(0, 1)

        self.ax1.set_ylim(-1.0, 1.1)
        self.ax1.set_xlim(0, 1)
        
        del self.pressure_mbn[:]
        del self.pressure_ruby[:]
        del self.step[:] 

        self.line1.set_data(self.step, self.pressure_mbn)
        self.line2.set_data(self.step, self.pressure_ruby) 

        return self.line1,self.line2,

    def SetupAnimation(self): 

        self.fig=plt.figure(figsize=(10,6))
        
        self.ax0 = self.fig.add_subplot() 
        self.ax1 = self.ax0.twinx()
    
        self.ax0.set_ylabel('Controller Output (bar)')

        self.ax1.set_ylabel('Ruby Pressure (GPa)')
        
        self.ax0.set_xlabel('Time (s)')
        self.ax0.grid()
        
        self.line1,= self.ax0.plot([],[], color='blue')
        self.line2,= self.ax1.plot([],[],color='red')
        
        self.line1.set_label('Controller Output')
        self.line2.set_label('Ruby Pressure')

        self.fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.3)
        
        self.fig.suptitle(f'PID Step Test \n Parameters: Pressure Set Point={self.PressureSetPoint} bar')  

        ani = animation.FuncAnimation(self.fig, self.run, self.generator, interval=1000, init_func=self.graph_init,
                              save_count=200000)
        
        ani.save('Animation.png')
        
        #plt.show()

    def run(self,data): 

        global time1
        t, y1, y2, y3= data
        
        time2=time.time() 
        experiment_time=round(time2-time1,3)  

        self.step.append(t)
        self.pressure_mbn.append(y1)
        self.pressure_ruby.append(y2) 
        self.mbn_real_pressure.append(y3)
        self.experiment_time.append(experiment_time)

        np.savetxt(f'StepTests/{self.FilePath}.csv',  [p for p in zip(self.step, self.pressure_mbn, self.pressure_ruby,self.mbn_real_pressure,self.experiment_time)],delimiter=',', fmt='%5.4f',
                header='{0:^5s},{1:^7s},{2:^9s},{3:^11s},{4:^13s}'.format('time(s)','Controler Setpoint (bar)','DAC Pressure (GPa)','Real Mbn Pressure (bar)','Experiment Time (s)',))
        
        print(f'-----------\nStep: {self.step[-1]}') 
        print(f'Controller Output: {self.pressure_mbn[-1]}')
        print(f'Mbn Real Pressure: {self.mbn_real_pressure[-1]}')
        print(f'Tempo de experimento: {experiment_time}') 
        print(f'Ruby Pressure: {self.pressure_ruby[-1]}\n-----------') 
        

        xmin1, xmax1 = self.ax0.get_xlim() 
        xmin2, xmax2 = self.ax1.get_xlim()

        ymin1, ymax1 = self.ax0.get_ylim() 
        ymin2, ymax2 = self.ax1.get_ylim()


        #Adjusting x limits
        if t >= xmax1:
            self.ax0.set_xlim(xmin1, 2*xmax1)
            self.ax0.figure.canvas.draw()
        
        if t >= xmax2:
            self.ax1.set_xlim(xmin2, 2*xmax2)
            self.ax1.figure.canvas.draw()


        #Adjusting y limits
        if y1 >= ymax1:
            self.ax0.set_ylim(ymin1, 1.5*ymax1)
            self.ax0.figure.canvas.draw()
        
        if y2 >= ymax2:
            self.ax1.set_ylim(ymin2, 1.5*ymax2)
            self.ax1.figure.canvas.draw()

        
        self.line1.set_data(self.step, self.pressure_mbn)
        self.line2.set_data(self.step, self.pressure_ruby) 

        return self.line1,self.line2,      

    def generator(self):  
        global i
        caput('EMA:B:PACE01:Control', 1)
        i=0

        self.SetPressure=self.CurrentMbnPressure 
        
        self.PressureRBVThread= DACPressureValueThread() 
        self.PressureRBVThread.start()
        self.PressureRBVThread.join()
        self.PressureRBV=self.PressureRBVThread.PressureValue  
        
        while True:
                         
            self.PressureRBVThread= DACPressureValueThread() 
            self.PressureRBVThread.start()
            self.PressureRBVThread.join()
            self.PressureRBV=self.PressureRBVThread.PressureValue  
             

            if i==10: 
                self.SetPressure=self.PressureSetPoint
                caput(pv_setpoint, self.SetPressure, wait=True) 
            
            self.RealMbnPressure=caget(pv_pressure_rbv)

            i+=1 
            yield i, self.SetPressure, self.PressureRBV, self.RealMbnPressure 


if __name__ == '__main__': 
    
    parser = argparse.ArgumentParser(description = 'Script to test the Step Response to the Pressure Controller.')

    #Positional Arguments
    parser.add_argument('pressure_set_point', type=float, help='The pressure (in bar) sent to the controller')
    parser.add_argument('current_mbn_pressure', type=float, help='The current membrane pressure')
    parser.add_argument('file_name', type=str, help='The name of the file where the data will be saved')
    args = parser.parse_args()

    PressureSetPoint=args.pressure_set_point
    CurrentMbnPressure=args.current_mbn_pressure
    FileName=args.file_name

    PIDClass(PressureSetPoint,CurrentMbnPressure,FileName)

