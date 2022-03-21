import os, sys
import time, glob

import random
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

import myokit


from simulator_myokit import Simulator


class Generator(Simulator):
    """
    
    """
    def __init__(self, model_path, protocol_def=None, pre_sim=2):
        super(Generator, self).__init__(model_path, protocol_def, pre_sim)
        '''
        if pre_sim==0 : No pre simulation
        if pre_sim==1 : pre simulation
        if pre_sim==2 : constant pacing pre simulation
        '''
 

    def gen_dataset(self, gen_params, datasetNo=1):
        '''
        type = 'AP' or 'I" 
        params = {
            'times': 1,                    
            'log_li' : [],
            'nData' : 10000,                         
            'dataset_dir' :   './dataset/,
            'data_file_name' :  'current,
            'scale' : 2,
        }  
        '''
        random.seed(datasetNo * 84)
        np.random.seed(datasetNo * 86)
    
        print("-----Dataset%d generation starts.-----"%(datasetNo))
        
        d = None              
        result_li = []
        param_li = []
        current_nData = 0
        
        simulation_error_count = 0
        with tqdm(total = gen_params['nData']) as pbar: 
            while (current_nData < gen_params['nData']):                
                g_fc_inv = np.random.uniform(0.0000001, 1, 8)                      
                if current_nData==0 and datasetNo==1:
                    g_fc_inv = np.ones(8)     
                g_fc = 1.0/g_fc_inv       # np.ones(6)
                
                Gfc = {                    
                    'ina.GNafc' : g_fc[0],
                    'inal.GNaLfc' : g_fc[1],
                    'ito.Gtofc' : g_fc[2],
                    'ical.PCafc' : g_fc[3],
                    'ikr.GKrfc' : g_fc[4],
                    'iks.GKsfc' : g_fc[5],
                    'ik1.GK1fc' : g_fc[6],
                    'if.Gffc' : g_fc[7]    
                } 
                self.set_simulation_params(Gfc)                
                log_li = ['membrane.V']
                if len(log_li)>0:
                    log_li = gen_params['log_li']
                try :
                    d = self.simulate( gen_params['times'], extra_log=gen_params['log_li'])                           
                    temp = []
                    for log in log_li :                                              
                        temp.append(d[log][::gen_params['scale']]) ###################
                    result_li.append( np.array(temp) ) 
                    param_li.append( g_fc_inv )
                    current_nData+=1                    
                except :
                    simulation_error_count += 1
    #                 print("There is a simulation error.")
                    continue
                pbar.update(1) 
  
        result_li = np.array(result_li, dtype=np.float32)        
        param_li = np.array(param_li, dtype=np.float32)        
                
        np.save(os.path.join(gen_params['dataset_dir'], 'times'), d['engine.time'].astype(np.float32)[::gen_params['scale']] )#######################3
        np.save(os.path.join(gen_params['dataset_dir'], 'voltage_protocol'), d['membrane.V'].astype(np.float32)[::gen_params['scale']] )#########################
        np.save(os.path.join(gen_params['dataset_dir'], "%s%d"%(gen_params['data_file_name'], datasetNo) ), result_li)
        np.save(os.path.join(gen_params['dataset_dir'], 'parameter%d'%(datasetNo) ), param_li )
        
        result_li = []
        param_li = []
            
        print("=====Dataset%d generation End.===== %d simulation errors occured.====="%(datasetNo, simulation_error_count))       
    
    
          
    

        
if __name__=='__main__':
    
    start_time = time.time()
   
   

    print("--- %s seconds ---"%(time.time()-start_time))
