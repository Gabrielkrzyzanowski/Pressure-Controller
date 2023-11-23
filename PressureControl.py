#!/usr/bin/env python3 

from epics import * 
import sys
import numpy as np
import spe_loader as sl
from threading import Thread
from PyQt5.QtCore import pyqtSignal, pyqtSlot, Qt, QThread, QObject,QThreadPool
from PyQt5.QtWidgets import QApplication, QMainWindow, QDialog,QWidget,QMessageBox
from PyQt5.QtGui import QFont
from PyQt5.QtGui import *
from PyQt5 import QtCore,uic 
import time 
import pyqtgraph as pg
from epics import caget,caput
import os

#Menu bar
from PyQt5.QtWidgets import QMenuBar, QMenu

#### import of my scripts 
from RubyPressure import CalculatePressure

### pvs controlador
file_path='/ibira/lnls/beamlines/ema/apps/LightField/Experiments_Files'
ruby_intesity=[]


prefix = 'EMA:B:PACE01:'
pv_controler_setpoint = prefix + 'Setpoint'
pv_controler_setpoint_rbv = pv_controler_setpoint + '_RBV'
pv_controler_slew = prefix + 'Slew'
pv_controler_pressure_rbv = prefix + 'Pressure_RBV'

pv_exposition_time= 'EMA:13LF1:cam1:AcquireTime'
pv_acquire_mode='EMA:13LF1:cam1:Acquire'
pv_file_number='EMA:13LF1:cam1:FileNumber'

pv_ROIX_Bin= 'EMA:13LF1:cam1:BinX'
pv_ROIY_Bin= 'EMA:13LF1:cam1:BinY' 
pv_ROIX_Size='EMA:13LF1:cam1:SizeX'
pv_ROIY_Size='EMA:13LF1:cam1:SizeY'

global i 
i=1

icon_green_led = './led-green-on.png'
icon_red_led = './led-red-on.png'



class GetData(QObject): ## Class that load the .SPE with the Ruby Spec Data
    
    params = pyqtSignal(list)
    stopped = pyqtSignal()
    
    def UpdateValues(self): 

        ruby_intesity=[]
        caput(pv_acquire_mode,1)
        
        
        experiment_files=os.listdir(file_path)
        experiment_files.sort()

        
        if len(experiment_files)<=3:

            while len(experiment_files)<3: 

                caput(pv_acquire_mode,1,wait=True)

                experiment_files=os.listdir(file_path)
                experiment_files.sort()

        ruby_data=sl.SpeFile(f'{file_path}/{experiment_files[1]}')          

        ruby_intesity.append(ruby_data.data[0][0][0]) 
        wave_lenght = ruby_data.wavelength

        RubySpec = ruby_intesity[0]
        WaveLenght = np.resize(wave_lenght, len(RubySpec))

        WaveLenght_ROI=WaveLenght[np.where((WaveLenght > 690) & (WaveLenght < 720))] 
        RubySpec_ROI= RubySpec[np.where((WaveLenght > 690) & (WaveLenght < 720))] 

        print(WaveLenght_ROI)
        os.remove(file_path +'/'+ experiment_files[0]) 

        if experiment_files[-1]=='test_970.spe': 
            caput(pv_file_number,1,wait=True) 
            caput(pv_acquire_mode,1,wait=True)
            time.sleep(0.2) 
            for f in os.listdir(file_path):
                if f=='test_001_spe': pass
                else: os.remove(os.path.join(file_path, f))

        
        self.params.emit([WaveLenght_ROI,RubySpec_ROI])
        self.stopped.emit()
        

class SetPressure(QObject): ##Class that uses my script to calculate the DAC pressure from the ruby spec
    
    pressure=pyqtSignal(float) 
    stopped = pyqtSignal()

    def __init__(self, *args, parent=None): 
        super(SetPressure,self).__init__(parent)   
        self.WaveLenght=args[0][0]
        self.RubySpec=args[0][1]

    def Setter(self):

        DacPressureValue = CalculatePressure(self.WaveLenght,self.RubySpec)
        DacPressureValue=round(DacPressureValue,3) 
        self.pressure.emit(DacPressureValue) 
        self.stopped.emit() 



class DACDisplay(QMainWindow): ##GUI Class 

    def __init__(self, *args, **kwargs):
        super(DACDisplay,self).__init__(*args, **kwargs) # Call the inherited classes __init__ method
        uic.loadUi('PressureControl.ui', self) # Load the ui file
        
        self.acquisition_frequency_value=1000
        self.WaveLenght=None
        self.RubySpec=None
        self.Acquisition=False
        self.Thread=QThread() 
        self._createMenuBar()
        self._createStatusBar() 
        self.LabelsSetup()
        self.MakeConnections() 
        self.DesignSetup()   
        self.PlotConfig() 

        self.RunSetMbnPressure()

        self.show()
    
    ##### UI Initial Config ##### 
    def _createMenuBar(self):
        menuBar = self.menuBar()
        # Creating menus using a QMenu object
        fileMenu = QMenu("&File", self)
        menuBar.addMenu(fileMenu)
        # Creating menus using a title
        editMenu = menuBar.addMenu("&Edit")
        helpMenu = menuBar.addMenu("&Help")

    def _createStatusBar(self):
        self.statusbar = self.statusBar()
        # Adding a temporary message
        self.statusbar.showMessage("Ready", 3000)

    def LabelsSetup(self): 
        self.currentExpositionTime=caget('EMA:13LF1:cam1:AcquireTime_RBV')
        self.exposition_time.setFont(QFont('Arial', 12))
        self.exposition_time.setText(f'{self.currentExpositionTime}')

        self.currentSlitWidth=caget('EMA:13LF1:cam1:LFEntranceWidth_RBV')
        self.Target_slit_width.setFont(QFont('Arial', 12))
        self.Target_slit_width.setText(f'{self.currentSlitWidth}')

        self.aquisition_frequency.setFont(QFont('Arial', 12))
        self.aquisition_frequency.setText(f'{self.acquisition_frequency_value}')

        self.currentGrating=caget('EMA:13LF1:cam1:LFGrating_RBV')
        self.Dac_list.setCurrentText(f'{self.currentGrating}') 

        self.pixmap_red= QPixmap(icon_red_led)
        self.pixmap_red= self.pixmap_red.scaledToWidth(30) 

        self.pixmap_green= QPixmap(icon_green_led)
        self.pixmap_green= self.pixmap_green.scaledToWidth(30) 

        #self.label_led.setPixmap(self.pixmap_green)

    def MakeConnections(self): 

        self.btn_set_grating.clicked.connect(self.SetGrating)
        self.btn_acquisition_start.clicked.connect(self.Timer)
        self.btn_acquisition_stop.clicked.connect(self.timer_killer)

        self.btn_set_dac.clicked.connect(self.ShowThreshold)
        self.btn_set_mbn_pressure.clicked.connect(self.SetNewMbnPressure)
        self.btn_exposition_time.clicked.connect(self.SetExpositionTime)
        self.btn_set_dac_pressure.clicked.connect(self.PIDRun)

        self.btn_set_slit_width.clicked.connect(self.SetSlitWidth)
        self.btn_aquisition_frequency.clicked.connect(self.ChangeTimer)

    def DesignSetup(self): 
        self.threshold_pressure.setStyleSheet('background-color: #fff4f4; border: 1px solid black;')
        self.threshold_pressure.setAlignment(QtCore.Qt.AlignCenter)
        self.threshold_pressure.setFont(QFont('Arial', 12))

        self.current_pressure.setStyleSheet('background-color: #fff4f4; border: 1px solid black;')
        self.current_pressure.setAlignment(QtCore.Qt.AlignCenter)
        self.current_pressure.setFont(QFont('Arial', 12)) 

        self.mbn_pressure.setStyleSheet('background-color: #9ab8ef; border: 1px solid black;')
        self.mbn_pressure.setAlignment(QtCore.Qt.AlignCenter)
        self.mbn_pressure.setFont(QFont('Arial', 12)) 

        self.mbn_pressure_variation.setStyleSheet('background-color: #9ab8ef; border: 1px solid black;')
        self.mbn_pressure_variation.setAlignment(QtCore.Qt.AlignCenter)
        self.mbn_pressure_variation.setFont(QFont('Arial', 12)) 

        self.btn_acquisition_start.setStyleSheet('background-color: #00ff7f')
        self.btn_acquisition_stop.setStyleSheet('background-color: #ff7474') 
        self.btn_aquisition_frequency.setStyleSheet('background-color: #00ff7f')
        
        
    def PlotConfig(self): 

        self.graphWidget = pg.PlotWidget() 
        self.graphWidget.addLegend()

        self.graphWidget.setBackground('w')

        pen = pg.mkPen(color=(255, 0, 0), width=2, style=QtCore.Qt.SolidLine) 
        self.pen2=pg.mkPen(color=(65,105,225), width=2, style=QtCore.Qt.DashDotLine)
        styles = {'color':'k', 'font-size':'14px'}

        self.graphWidget.setLabel('left', 'Counts', **styles)
        self.graphWidget.setLabel('bottom', 'Wave Lenght (nm)', **styles) 

        self.style = pg.PlotDataItem(pen=None, symbol='o')  
        self.graphWidget.plotItem.legend.addItem(self.style, 'label') 

        self.graphWidget.showGrid(x=True, y=True)

        self.pseudo_line=self.graphWidget.plot([], [], pen=None, symbol='o')
        self.data_line=self.graphWidget.plot([], [], pen=pen, name='Ruby Spec') 
        
        self.Graph_widget.addWidget(self.graphWidget)

    ######### Get Data ######### 
   
    def Timer(self):
        self.Acquisition=True         
        self.timer = QtCore.QTimer()
        self.timer.setInterval(self.acquisition_frequency_value)
        self.timer.timeout.connect(self.DataLoader)
        self.timer.timeout.connect(self.PressureLoader)
        self.timer.start() 

    def timer_killer(self):
        while caget('EMA:13LF1:cam1:AcquireBusy')!=0: 
            pass
        caput(pv_acquire_mode,0,wait=True)
        self.timer.stop()
        self.Acquisition=False

    def ChangeTimer(self):
        if self.Acquisition == False:
             
            self.acquisition_frequency_value_backup=self.acquisition_frequency_value
            self.acquisition_frequency_value=int(self.aquisition_frequency.text())

            if self.acquisition_frequency_value < 2*self.currentExpositionTime:
                self.show_warning_acquisition_frequency()
                self.acquisition_frequency_value=self.acquisition_frequency_value_backup
                self.LabelsSetup() 
        elif self.Acquisition == True:
            self.show_warning_changes()


    def DataLoader(self): 
        if self.Acquisition==True:
            
            self.DataObject= GetData()
            self.DataObject.moveToThread(self.Thread) 
            self.Thread.started.connect(self.DataObject.UpdateValues)

            self.DataObject.params.connect(self.UpdateData)
            self.DataObject.stopped.connect(self.Thread.quit) 
        
            self.Thread.start()


    @pyqtSlot(list)
    def UpdateData(self,params): 

        self.WaveLenght=params[0] 
        self.RubySpec=params[1]  
        print(self.WaveLenght,self.RubySpec)
        self.graphWidget.plotItem.legend.removeItem(self.style)
        self.graphWidget.plotItem.legend.addItem(self.style,
                                                f'Peak position: {self.WaveLenght[np.where(self.RubySpec==max(self.RubySpec))][0]:.2f} nm')


        if self.y_max.text() and self.y_min.text():
            self.graphWidget.setYRange(int(self.y_min.text()),int(self.y_max.text()))
        if self.x_max.text() and self.x_min.text(): 
            self.graphWidget.setXRange(int(self.x_min.text()),int(self.x_max.text()))
        
        self.pseudo_line.setData(x=[self.WaveLenght[np.where(self.RubySpec==max(self.RubySpec))][0]] , y=[np.max(self.RubySpec)]) 

        self.data_line.setData(self.WaveLenght, self.RubySpec)
        self.SaturationCorrection(self.WaveLenght,self.RubySpec)
        self.FindSaturatedPoints(self.WaveLenght,self.RubySpec)
    

    def PressureLoader(self): 
        
        if self.WaveLenght is not None: 
    
            self.PressureObject=SetPressure((self.WaveLenght,self.RubySpec))
            self.PressureObject.moveToThread(self.Thread)
            self.Thread.started.connect(self.PressureObject.Setter) 
            
            self.PressureObject.pressure.connect(self.PrintPressure) 
            self.PressureObject.stopped.connect(self.Thread.quit)
            self.Thread.start() 

    ######### Display Information ######### 

    def RunSetMbnPressure(self): 
        self.mbntimer = QtCore.QTimer()
        self.mbntimer.setInterval(200)
        self.mbntimer.timeout.connect(self.SetMbnPressure)
        self.mbntimer.start() 


    def SetMbnPressure(self): 

        self.current_mbn_pressure= caget(pv_controler_pressure_rbv) 
        self.current_mbn_pressure_variation=caget(pv_controler_slew)
        self.mbn_pressure.setText(f'{self.current_mbn_pressure:.3f} bar')
        self.mbn_pressure_variation.setText(f'{self.current_mbn_pressure_variation} bar/s')

  
    @pyqtSlot(float)
    def PrintPressure(self,DacPressureValue):
        global PressureDac 
        PressureDac=DacPressureValue 

        self.DacPressureValue=DacPressureValue
        self.current_pressure.setText(f'{self.DacPressureValue:.3f} GPa') 

    
    def SaturationCorrection(self,WaveLenght,RubySpec): 
        Max_intensity= np.max(RubySpec) 
        Max_intensity_index= np.where(RubySpec==Max_intensity)[0][0] 

        WaveLenghtAdjusted=WaveLenght[Max_intensity_index -45 :Max_intensity_index+20]
        RubyIntensityAdjusted=RubySpec[Max_intensity_index -45 : Max_intensity_index+20] 

        grady=np.gradient(RubyIntensityAdjusted) 
        gradx=np.gradient(WaveLenghtAdjusted) 

        dydx=grady/gradx 

        threshold = 0
        threshold_array = np.full(len(dydx), fill_value = threshold)

        if True in np.isclose(dydx, threshold_array, rtol = 1e-8): 
            print(WaveLenght[np.where(np.isclose(dydx, threshold_array, rtol = 1e-8) == True)])
            


    def SetGrating(self): 

        if self.Acquisition == False:
             
            if self.Grating_list.currentText()=='Density: 2400 g/mm ; Holographic: VIS': 
                caput('EMA:13LF1:cam1:LFGrating', '[h-vis,2400][0][0]')

            elif self.Grating_list.currentText()=='Density: 1800 g/mm ; Blaze: 500 nm':  
                caput('EMA:13LF1:cam1:LFGrating', '[500nm,1800][1][0]')

            elif self.Grating_list.currentText()=='Density: 1200 g/mm ; Blaze: 500 nm': 
                caput('EMA:13LF1:cam1:LFGrating', '[500nm,1200][2][0]')
        
            else:
                self.show_warning_grating()

        elif self.Acquisition == True:
            self.show_warning_changes()


    def ShowThreshold(self): 

        if self.Dac_list.currentText()=='350d001': 
            self.threshold_pressure.setText(f'{18.0} bar')

        elif self.Dac_list.currentText()=='350d003': 
            self.threshold_pressure.setText(f'{15.0} bar')

        elif self.Dac_list.currentText()=='250d004': 
            self.threshold_pressure.setText(f'{15.0} bar')

        elif self.Dac_list.currentText()=='630d00#': 
            self.threshold_pressure.setText(f'{15.0} bar')

        else: 
            self.show_warning_dac()

    def show_warning_changes(self,wtitle="Acquisition Running",
                        txt="You must stop your acquisition before any change.",infotxt = None): ## Show warning when no dac is selected 
        
        wrn = QMessageBox()
        wrn.setWindowTitle(wtitle)
        wrn.setIcon(QMessageBox.Warning)
        wrn.setText(txt)
        wrn.setInformativeText(infotxt)
        wrn.setStandardButtons(QMessageBox.Ok)
        
        retval = wrn.exec_()

        return retval  


    def show_warning_grating(self,wtitle="No Grating selected",
                        txt="Select the Grating to proceed.",infotxt = None): ## Show warning when no dac is selected 
        
        wrn = QMessageBox()
        wrn.setWindowTitle(wtitle)
        wrn.setIcon(QMessageBox.Warning)
        wrn.setText(txt)
        wrn.setInformativeText(infotxt)
        wrn.setStandardButtons(QMessageBox.Ok)
        
        retval = wrn.exec_()

        return retval  


    def show_warning_dac(self,wtitle="No Dac selected",
                        txt="Select the dac to proceed.",infotxt = None): ## Show warning when no dac is selected 
        
        wrn = QMessageBox()
        wrn.setWindowTitle(wtitle)
        wrn.setIcon(QMessageBox.Warning)
        wrn.setText(txt)
        wrn.setInformativeText(infotxt)
        wrn.setStandardButtons(QMessageBox.Ok)
        
        retval = wrn.exec_()

        return retval 

    def show_warning_acquisition_frequency(self,wtitle="Acquisition Frequency time too low",
                        txt="You must set an Acquisition Frequency value at least\ntwo times bigger than the current Exposition Time..",infotxt = None): ## Show warning when no dac is selected 
        
        wrn = QMessageBox()
        wrn.setWindowTitle(wtitle)
        wrn.setIcon(QMessageBox.Warning)
        wrn.setText(txt)
        wrn.setInformativeText(infotxt)
        wrn.setStandardButtons(QMessageBox.Ok)
        
        retval = wrn.exec_()

        return retval 
  

    def show_warning_mbn_pressure(self,wtitle="No Target Pressure",
                        txt="Set a value to the Target Pressure to proceed.",infotxt = None): ## Show warning when no dac is selected 
        
        wrn = QMessageBox()
        wrn.setWindowTitle(wtitle)
        wrn.setIcon(QMessageBox.Warning)
        wrn.setText(txt)
        wrn.setInformativeText(infotxt)
        wrn.setStandardButtons(QMessageBox.Ok)
        
        retval = wrn.exec_()

        return retval 


    def show_warning_mbn_pressure2(self,wtitle="Target Pressure too low",
                        txt="The Target Pressure value is lower than the Current Pressure \n               Set a higher pressure value to proceed",
                        infotxt = None): 
        
        wrn = QMessageBox()
        wrn.setWindowTitle(wtitle)
        wrn.setIcon(QMessageBox.Warning)
        wrn.setText(txt)
        wrn.setInformativeText(infotxt)
        wrn.setStandardButtons(QMessageBox.Ok)
        
        retval = wrn.exec_()

        return retval 
    

    def show_warning_exposition_time(self,wtitle="No Exposition Time",
                        txt="Set a value to the Exposition Time to proceed.",infotxt = None): ## Show warning when no dac is selected 
        
        wrn = QMessageBox()
        wrn.setWindowTitle(wtitle)
        wrn.setIcon(QMessageBox.Warning)
        wrn.setText(txt)
        wrn.setInformativeText(infotxt)
        wrn.setStandardButtons(QMessageBox.Ok)
        
        retval = wrn.exec_()

        return retval 


    def show_warning_slit(self,wtitle="Slit Width value out of bounds",
                        txt="Set an aceptable value to the Slit Width to proceed. \nRange: 10 um to 3000 um",infotxt = None): ## Show warning when no dac is selected 
        
        wrn = QMessageBox()
        wrn.setWindowTitle(wtitle)
        wrn.setIcon(QMessageBox.Warning)
        wrn.setText(txt)
        wrn.setInformativeText(infotxt)
        wrn.setStandardButtons(QMessageBox.Ok)
        
        retval = wrn.exec_()

        return retval 

    
    ######### Active Processes and Calculations #########   

    def PIDRun(self): 
        
        self.ThreadPID=QThread()

        self.ShowThreshold() #In case there is no DAC selected 
        
        self.Pressure_alarm=int(self.Pressure_Alarm.text())
        self.Pressure_setpoint=int(self.Target_Pressure.text())

        self.PIDObject= PIDClass(self.Pressure_setpoint, self.Pressure_alarm, FilePath='Tests/rubytest1.spe')
        self.PIDObject.moveToThread(self.ThreadPID) 
        self.ThreadPID.started.connect(self.PIDObject.generator)

        self.PIDObject.stopped.connect(self.ThreadPID.quit) 
        
        self.ThreadPID.start()


    def SetNewMbnPressure(self):  
        
        if self.set_mbn_pressure.text() == '': 
            self.show_warning_mbn_pressure() 
        elif float(self.set_mbn_pressure.text()) <= self.current_mbn_pressure : 
            self.show_warning_mbn_pressure2()
        else:
            caput('EMA:B:PACE01:Control', 1)
            print(f'Setting Membrane Pressure to: {self.set_mbn_pressure.text()} bar')
            caput(pv_controler_setpoint, self.set_mbn_pressure.text(), wait=True)
            
    
    def SetExpositionTime(self): 
        
        if self.exposition_time.text() == '': 
            self.show_warning_exposition_time()
        else:
            self.exposition_time_value=float(self.exposition_time.text())
            if self.acquisition_frequency_value < 2*self.exposition_time_value:
                self.show_warning_acquisition_frequency()
                self.LabelsSetup() 
            else:
                caput(pv_exposition_time, f'{self.exposition_time_value/1000}')
                self.currentExpositionTime = self.exposition_time_value


    def SetSlitWidth(self): 
        
        if self.Acquisition == False:
             
            if self.Target_slit_width.text() == '': 
                self.show_warning_slit() 
            elif int(self.Target_slit_width.text()) < 10 : 
                self.show_warning_slit()
            elif int(self.Target_slit_width.text()) > 3000 : 
                self.show_warning_slit()
            else: 
                caput('EMA:13LF1:cam1:LFEntranceWidth', self.Target_slit_width.text()) 

        elif self.Acquisition == True:
            self.show_warning_changes()



    def FindSaturatedPoints(self, WaveLenght,RubySpec): 
        
        Max_intensity= np.max(RubySpec) 
        Max_intensity_index= np.where(RubySpec==Max_intensity)[0][0]

        Resized_WaveLenght=WaveLenght[Max_intensity_index -20 :Max_intensity_index+20]
        Resized_RubySpec=RubySpec[Max_intensity_index -20 : Max_intensity_index+20]

        grady=np.gradient(Resized_RubySpec) 
        gradx=np.gradient(Resized_WaveLenght)

        dydx=grady/gradx

        threshold = 0
        thresh_array = np.full(len(dydx), fill_value = threshold) 

        SaturatedPoints= WaveLenght[np.where(np.isclose(dydx, thresh_array, rtol = 1e-9) == True)] 

        self.Saturation_indicator.connectionStateChanged(1)  if len(SaturatedPoints) >= 2 and Max_intensity_index in SaturatedPoints  else 0 


class PIDClass(QObject):

    stopped = pyqtSignal()

    def __init__(
            self,
            PressureSetPoint: float = None,
            Kp: float = None,
            Kc: float = None, 
            Kd: float= None, 
            PressureAlarm: float = None,
            FilePath: str = None
            ) -> None: 
        
        super(PIDClass,self).__init__()
    
        self.PressureSetPoint=PressureSetPoint

        self.Kp=Kp if Kp!=None else 0.01 #Calibrated Value 
        self.Kc=Kc if Kc!=None else 0.01 #Calibrated Value 
        self.Kd=Kd if Kd!=None else 0.01 #Calibrated Value 

        self.PressureAlarm=PressureAlarm
        self.FilePath=FilePath 


    def generator(self):  
        
        i=0 

        if 'PressureDac' not in globals() : 
            print('NO PRESSURE VALUE FOUND!!!!!! \nClosing the thread...')
            self.stopped.emit()
        else: 

            self.PressureRBV=PressureDac

            controler = self.LoopPID(self.Kp,self.Kc,self.Kd) 
            controler.send(None)
            self.SetPressure = controler.send((i,self.PressureRBV, self.PressureSetPoint))
        
            while abs(self.PressureSetPoint - self.PressureRBV) >= 1e-6:
                
                self.PressureRBV = PressureDac

                controler = self.LoopPID(self.Kp,self.Kc,self.Kd) 
                controler.send(None)
                
                self.SetPressure = controler.send((i,self.PressureRBV, self.PressureSetPoint))
                print(f'SetPressure: {self.SetPressure}')

                if self.SetPressure >= self.PressureAlarm: 
                    raise SystemExit('ALARM PRESSURE EXCEED!!!!!! \nClosing the script...')
                else:   
                    caput('pv_setpoint', self.SetPressure, wait=True)
    
                while (caget('pv_pressure_rbv') - self.SetPressure) >= 1e-3:
                    pass
                  
                i+=1 

                time.sleep(1) 

            self.stopped.emit() 
        

    def LoopPID(self,Kp,Ki,Kd,Pv_threshold=0): #Here Pv_threshold is the inicial pressure in the membrane
        #Initializing
        e_prev=0 
        t_prev=-100
        SetPressure=0
        P,I,D=0,0,0

        while True:
            # yield SetPressure, wait for new t, PressureRBV and PressureSetPoint
            t , PressureRBV , PressureSetPoint = yield SetPressure
            
            error= PressureSetPoint - PressureRBV #Error 

            P += Kp*error
            I += Ki*error*(t-t_prev)
            D += Kd*(error-e_prev)/(t-t_prev)

            SetPressure=Pv_threshold+P+I+D

            #update stored data for next iteration
            e_prev=error 
            t_prev=t 


if __name__ == "__main__": 

    app = QApplication(sys.argv) 
    windows = DACDisplay() # Create an instance of our class
    sys.exit(app.exec()) #Start the apllication 
