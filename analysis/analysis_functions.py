# -*- coding: utf-8 -*-
"""

@author: Karl Mayer
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Analysis functions for logical T benchmarking

def analyze_T_results(job_data, postselect=False):
    
    # read in sequence length and sequence reps
    #seq_len = list(set([job.program.parameters[0] for job in batch.jobs]))
    seq_len = list(set([job_data[key]['parameters'][0] for key in job_data]))
    seq_len.sort()
    
    succ_probs = {L:[] for L in seq_len}
    succ_stds = {L:[] for L in seq_len}
    ps_shots_dict = {L:[] for L in seq_len}
    
    for key in job_data:
        if 'results' in job_data[key]['results']:
            params = job_data[key]['parameters']
            L = params[0]
            jr = job_data[key]['results']['results']
            outcomes = {'0':0, '1':0}
            for i in range(len(jr['c_log'])):
                # postselect on init
                if jr['init'][i] == '0':
                    
                    if postselect == False:
                        b_str = jr['c_log'][i]
                        if b_str in outcomes:
                            outcomes[b_str] += 1
                        else:
                            outcomes[b_str] = 1
                    
                    # postselect on T gadget syndromes
                    elif postselect == True:
                        syn = ''
                        for t in range(1,L+1):
                            syn = syn + jr[f'syn_measT{t}'][i]
                        if syn.count('1') == 0:
                            b_str = jr['c_log'][i]
                            if b_str in outcomes:
                                outcomes[b_str] += 1
                            else:
                                outcomes[b_str] = 1
                    
            try:
                exp_out = params[2]
            except:
                exp_out = '0'
                
            ps_shots = sum(outcomes.values())
            ps_shots_dict[L].append(ps_shots)
            
            p = outcomes[exp_out]/ps_shots
            p_std = np.sqrt(p*(1-p)/ps_shots)
            succ_probs[L].append(p)
            succ_stds[L].append(p_std)
    
    data = {'succ':succ_probs, 'stds':succ_stds}
    if postselect == True:
        data['ps_shots_dict'] = ps_shots_dict
    
    return data
    
    
def plot_T_data(data, method=2, save=False,ylim=None,
             colors=['blue', 'blue', 'green', 'green'],
             labels=['', '', '', ''],
             linestyles=['-','--','-','--'],
             formats=['o','d','o','d'],
             filename='T_bench_plot_new.pdf'):
    """ data: dict or list of data dicts """
    
    def fit_func(x, a, f):
        return a*f**x + 1/2
    
    if type(data) == dict:
        data = [data]
        
    if type(method) == int:
        method = [method]
        
    y_list = []
    yerr_list = []
    popt_list = []
    perr_list = []
    yfit_list = []
    
    
    x = list(data[0]['succ'].keys())
    xfit = np.linspace(x[0], x[-1], 100)
    for D in data:

        y = [np.mean(D['succ'][L]) for L in x]
        shot_var = [sum([s**2 for s in D['stds'][L]])/len(D['stds'][L])**2 for L in x]
        circ_var = [np.std(D['succ'][L])**2/len(D['succ'][L]) for L in x]
        yerr = [np.sqrt(shot_var[j] + circ_var[j]) for j in range(len(x))]
        
        y_list.append(y)
        yerr_list.append(yerr)

        popt, pcov = curve_fit(fit_func, x, y)
        perr = np.sqrt(np.diag(pcov))
        yfit = fit_func(xfit, *popt)
        
        popt_list.append(popt)
        perr_list.append(perr)
        yfit_list.append(yfit)
    
    for i, y in enumerate(y_list):
        plt.errorbar(x, y, yerr=yerr_list[i], fmt=formats[i], color=colors[i])
        plt.plot(xfit, yfit_list[i], linestyles[i], color=colors[i], label=labels[i])
    plt.xlabel('Sequence Length (number of T gates)', fontsize=12)
    plt.xticks(ticks=x, labels=x, fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylabel('Survival Probability', fontsize=12)
    if save == False:
        plt.title('Logical non-FT T gate, MeasDecode')   
    
    #plt.legend(fontsize=12)
    plt.legend(fontsize=12)    
    plt.ylim(ylim)
    if save == True:
        plt.savefig(filename, format='pdf')
    plt.show()
    
    for i, popt in enumerate(popt_list):
        f_est, f_std = popt[-1], perr_list[i][-1]
        log_err, log_err_std = (1-f_est)/2, f_std/2
        F_avg, F_avg_std = (2*(1-log_err)+1)/3, 2*log_err_std/3

        print(f'T method {method[i]}')
        print(f'Prob of logical Z error: {round(log_err,4)} +/- {round(log_err_std,4)}')
        print(f'Average Fidelity:        {round(F_avg,4)} +/- {round(F_avg_std,4)}\n')


# Analysis functions for logical QFT 



def analyze_QFT_results(job_data, method='ancilla-assisted', postselect=False):
    
    data = {'comp':{'succ':{}, 'stds':{}}, 'fourier':{'succ':{}, 'stds':{}}}
    
    if method == 'ancilla-assisted':
        init_str = '0000'
    elif method == 'recursive-teleportation':
        init_str = '000'
    
    for key in job_data:
        if 'results' in job_data[key]['results']:
            params = job_data[key]['parameters']
            jr = job_data[key]['results']['results']
            outcomes = {}
            for i in range(len(jr['c'])):
                
                # postselect on init string
                if jr['init'][i] == init_str:
                    
                    # postselect on T gadget measurements
                    if postselect == False:
                        b_str = jr['c'][i]
                        if b_str in outcomes:
                            outcomes[b_str] += 1
                        else:
                            outcomes[b_str] = 1
                    
                    elif postselect == True:
                        syn = ''
                        for t in range(1,14):
                            syn = syn + jr[f'syn_measT{t}'][i]
                            
                        if syn.count('1') == 0:
                            b_str = jr['c'][i]
                            if b_str in outcomes:
                                outcomes[b_str] += 1
                            else:
                                outcomes[b_str] = 1

            # success probability
            exp_out = params[1]
            state = exp_out
            ps_shots = sum(outcomes.values())
            p = outcomes[exp_out]/ps_shots
            p_std = np.sqrt(p*(1-p)/ps_shots)
            
            data[params[0]]['succ'][state] = p
            data[params[0]]['stds'][state] = p_std
            
                
    return data


def plot_QFT_data(data, machine='H2-1', method='ancilla-assisted'):
    
    x = ['000', '001', '010', '011', '100', '101', '110', '111']
    x1, x2 = [], []
    y1, y2, yerr1, yerr2 = [], [], [], []
    for state in x:
        if state in data['comp']['succ']:
            x1.append(state)
            y1.append(data['comp']['succ'][state])
        if state in data['fourier']['succ']:
            x2.append(state)
            y2.append(data['fourier']['succ'][state])
        if state in data['comp']['stds']:
            yerr1.append(data['comp']['stds'][state])
        if state in data['fourier']['stds']:
            yerr2.append(data['fourier']['stds'][state])
    
    if y1 != []:
        f1 = np.mean(y1)
        f_std1 = np.sqrt(sum([s**2 for s in yerr1]))/len(y1)
    if y2 != []:
        f2 = np.mean(y2)
        f_std2 = np.sqrt(sum([s**2 for s in yerr2]))/len(y2)


    #w = 0.4
    if y1 != []:
        plt.bar(x1, y1, yerr=yerr1, label=machine)
        plt.xlabel('Input State', fontsize=12)
        plt.ylabel('State Fidelity', fontsize=12)
        plt.title(f'N=3 {method} QFT Comp Basis')
        plt.ylim(0,1)
        plt.legend()
        #plt.savefig('QFT1_comp.pdf', format='pdf')
        plt.show()
    
    if y2 != []:
        plt.bar(x2, y2, yerr=yerr2, label=machine)
        plt.xlabel('Input State', fontsize=12)
        plt.ylabel('State Fidelity', fontsize=12)
        plt.title(f'N=3 {method} QFT Fourier Basis')
        plt.ylim(0,1)
        plt.legend()
        #plt.savefig('QFT1_four.pdf', format='pdf')
        plt.show()


    if y1 != []:
        print(f'{method}, comp basis: {round(f1,3)} +/- {round(f_std1, 3)}')
    if y2 != []:
        print(f'{method}, four basis: {round(f2,3)} +/- {round(f_std2, 3)}')
        
    f_lo = (8*(f1+f2-1)+1)/9
    f_lo_std = 8*np.sqrt(f_std1**2+f_std2**2)/9
    
    print(f'{method}, F_avg_lo: {round(f_lo,3)} +/- {round(f_lo_std, 3)}')
        
        
def plot_QFT_comparison_data(data1, data2,
                         labels=['', ''],
                         colors=['blue','green'],
                         method='ancilla-assisted',
                         save=False):
    
    xlabels = ['000', '001', '010', '011', '100', '101', '110', '111']
    y1c, y2c, y1f, y2f = [], [], [], []
    y1c_err, y2c_err, y1f_err, y2f_err = [], [], [], []
    for state in xlabels:
        y1c.append(data1['comp']['succ'][state])
        y2c.append(data2['comp']['succ'][state])
        y1f.append(data1['fourier']['succ'][state])
        y2f.append(data2['fourier']['succ'][state])
        y1c_err.append(data1['comp']['stds'][state])
        y2c_err.append(data2['comp']['stds'][state])
        y1f_err.append(data1['fourier']['stds'][state])
        y2f_err.append(data2['fourier']['stds'][state])


    w = 0.4
    
    x = np.array(list(range(8)))
    plt.bar(x-w/2, y1c, width=w, yerr=y1c_err, label=f'{labels[0]}', color=colors[0])
    plt.bar(x+w/2, y2c, width=w, yerr=y2c_err, label=f'{labels[1]}', color=colors[1])
    plt.xticks(ticks=x, labels=xlabels)
    plt.xlabel('Input State', fontsize=12)
    plt.ylabel('State Fidelity', fontsize=12)
    if save == False:
        plt.title(f'N=3 QFT {method} Comp Basis')
    plt.ylim(0,1.1)
    plt.legend(loc=1)
    if save == True:
        plt.savefig(f'QFT_{method}_comp_new.pdf', format='pdf')
    plt.show()
    
    x = np.array(list(range(8)))
    plt.bar(x-w/2, y1f, width=w, yerr=y1f_err, label=f'{labels[0]}', color=colors[0])
    plt.bar(x+w/2, y2f, width=w, yerr=y2f_err, label=f'{labels[1]}', color=colors[1])
    plt.xticks(ticks=x, labels=xlabels)
    plt.xlabel('Input State', fontsize=12)
    plt.ylabel('State Fidelity', fontsize=12)
    if save == False:
        plt.title(f'N=3 QFT {method} Fourier Basis')
    plt.ylim(0,1.1)
    plt.legend(loc=1)
    if save == True:
        plt.savefig(f'QFT_{method}_four_new.pdf', format='pdf')
    plt.show()