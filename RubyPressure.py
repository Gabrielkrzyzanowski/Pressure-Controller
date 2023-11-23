#!/usr/bin/env python3 

"""
A Python script that calculates the value of the pressure in a Diamond Anvil Cell (DAC). 

The script get as an input a difraction spectrum from a Ruby inside the DAC, make a fit of this spectrum and use the formula 
derived for Shen et.all (https://doi.org/10.1080/08957959.2020.1791107) to relate the shift in the spectrum to the actual pressure
in the cell. 

Parameters
---------- 
wave_lenght: Array
    Array containing the wave lenghts in the spectrum of the ruby 
ruby_intensity: Array 
    Array containing the intensity plot in function of the wave lenght of the ruby spectrum

Ouput
---------- 
pressure: float
    The real value of the pressure in the DAC

"""
import numpy as np 
from lmfit.models import LinearModel, PseudoVoigtModel 
import random
from lmfit import models
from scipy.signal import find_peaks 


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

    peak_parameters={'prominence':1, 'width':3, 'height':4000}
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

