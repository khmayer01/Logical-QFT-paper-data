# -*- coding: utf-8 -*-
"""

@author: Karl Mayer
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# Analysis functions for logical TQRB

def analyze_TQRB_results(job_data):
    
    surv_prob = {}
    ps_shots = {}
    for key in job_data:
        params = job_data[key]['parameters']
        L = params['seq_len']
        if 'surv_state' in params:
            surv_state = params['surv_state']
        else:
            surv_state = '00'
            
        if L not in surv_prob:
            surv_prob[L] = []
            ps_shots[L] = []

        # compute outcomes
        outcomes = {}
        if 'results' in job_data[key]['results']:
            raw_results = job_data[key]['results']['results']['c']
            for b_str in ['00', '01', '10', '11']:
                if b_str in raw_results:
                    outcomes[b_str] = raw_results.count(b_str)

        # compute survival probability
        if outcomes != {}:
            shots = sum(outcomes.values())
            ps_shots[L].append(shots)
            if surv_state in outcomes:
                p = outcomes[surv_state]/shots
            else:
                p = 0

            surv_prob[L].append(p)
            
    # compute survival standard deviations
    surv_stds = {2:[], 6:[], 10:[], 14:[]}
    for L in surv_stds:
        for i in range(10):
            # rule of 2 for error bar for point near 100%
            p_adj = 100*surv_prob[L][i]/102
            std = np.sqrt(p_adj*(1-p_adj)/102)
            surv_stds[L].append(std)    
                
    data = {'surv_prob':surv_prob, 'ps_shots':ps_shots, 'surv_stds':surv_stds}
        
    return data


def plot_TQRB_data(data1, data2, save=False,
                   labels=['H1-1', 'H2-1'],
                   colors=['b', 'g'],
                   ylim=(0.85, 1.01)):

    
    def fit_func(x, a, f):
        return a*f**x + 1/4

    seq_len = [2, 6, 10, 14]
    
    x = seq_len
    avg_surv_probs1 = {L:np.mean(data1['surv_prob'][L]) for L in x}
    avg_surv_probs2 = {L:np.mean(data2['surv_prob'][L]) for L in x}
    surv_stds1 = data1['surv_stds']
    #ps_shots1 = data1['ps_shots']
    ps_shots2 = data2['ps_shots']
    y1 = [avg_surv_probs1[L] for L in x]
    y2 = [avg_surv_probs2[L] for L in x]
    
    shot_var1 = [sum([s**2 for s in surv_stds1[L]])/100 for L in x]
    #shot_var1 = [avg_surv_probs1[L]*(1-avg_surv_probs1[L])/sum(ps_shots1[L]) for L in x]
    # compute error bar using "rule of 2" for L=2 point
    #p_adj = (999/1002)
    #shot_var1[0] = (p_adj*(1-p_adj)/1000)
    shot_var2 = [avg_surv_probs2[L]*(1-avg_surv_probs2[L])/sum(ps_shots2[L]) for L in x]
    
    
    circ_var1 = [np.std([data1['surv_prob'][L][s] for s in range(10)])**2/10 for L in x]
    circ_var2 = [np.std([data2['surv_prob'][L][s] for s in range(10)])**2/10 for L in x]
    
    
    yerr1 = [np.sqrt(shot_var1[j]+circ_var1[j]) for j in range(4)]
    yerr2 = [np.sqrt(shot_var2[j]+circ_var2[j]) for j in range(4)]
    
    
    p0 = [0.75, 0.99]
    bounds = ([0, 0],[1, 1])
    
    popt1, pcov1 = curve_fit(fit_func, x, y1, p0=p0, bounds=bounds)
    popt2, pcov2 = curve_fit(fit_func, x, y2, p0=p0, bounds=bounds)
    
    xfit = np.linspace(x[0], x[-1], 100)
    yfit1 = fit_func(xfit, *popt1)
    yfit2 = fit_func(xfit, *popt2)
    
    plt.errorbar(x, y1, yerr=yerr1, fmt='o', color=colors[0])
    plt.errorbar(x, y2, yerr=yerr2, fmt='o', color=colors[1])
    
    plt.plot(xfit, yfit1, '-', color=colors[0], label=labels[0])
    plt.plot(xfit, yfit2, '-', color=colors[1], label=labels[1])
    
    plt.xlabel('Sequence Length (number of Cliffords)', fontsize=12)
    plt.xticks(ticks=x, labels=x, fontsize=12)
    plt.ylabel('Survival Probability', fontsize=12)
    plt.legend(fontsize=12)
    plt.ylim(ylim)
    if save == True:
        plt.savefig('TQRB_plot.pdf', format='pdf')
    plt.show()
    
    f1, f_std1 = popt1[1], np.sqrt(np.diag(pcov1))[1]
    f2, f_std2 = popt2[1], np.sqrt(np.diag(pcov2))[1]
    
    # 1.78073 is the average TQ count per Clifford in the construction used
    F_avg1 = 1 - 3*((1-f1)/1.78073)/4
    F_avg2 = 1 - 3*((1-f2)/1.78073)/4
    
    F_avg_std1 = (3/4)*f_std1/1.78073
    F_avg_std2 = (3/4)*f_std2/1.78073
    
    
    print(f'H1-1 Average Fidelity = {round(F_avg1, 4)} +/- {round(F_avg_std1, 4)}')
    print(f'H2-1 Average Fidelity = {round(F_avg2, 4)} +/- {round(F_avg_std2, 4)}')

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
    
    
# Analysis functions for logical control-T

def analyze_CT_results(job_data, postselect=False):

    data = {}

    for key in job_data:
        if 'results' in job_data[key]['results']:
            params = job_data[key]['parameters']
            jr = job_data[key]['results']['results']
            outcomes = {}
            for i in range(len(jr['c'])):

                # postselect on init
                if jr['init'][i] == '000':

                    if postselect == True:
                        if [jr[f'syn_measT{j}'][i] for j in range(1,7)] == ['000' for j in range(6)]:
                            b_str = jr['c'][i]
                            if b_str in outcomes:
                                outcomes[b_str] += 1
                            else:
                                outcomes[b_str] = 1

                    elif postselect == False:
                        b_str = jr['c'][i]
                        if b_str in outcomes:
                            outcomes[b_str] += 1
                        else:
                            outcomes[b_str] = 1

            data[params] = outcomes
    
    return data


def plot_CT_data(data, save=False):
    
    F1s = [data[state]['00']/sum(data[state].values()) for state in ['+0', '+1', '-0', '-1']]
    F1_shots = [sum(data[state].values()) for state in ['+0', '+1', '-0', '-1']]
    F1stds = [np.sqrt(p*(1-p)/F1_shots[i]) for i, p in enumerate(F1s)]

    F2s = [data[state]['00']/sum(data[state].values()) for state in ['0+', '0-', '1+', '1-']]
    F2_shots = [sum(data[state].values()) for state in ['0+', '0-', '1+', '1-']]
    F2stds = [np.sqrt(p*(1-p)/F2_shots[i]) for i, p in enumerate(F2s)]


    F1 = np.mean(F1s)
    F2 = np.mean(F2s)

    F1std = np.sqrt(sum([s**2 for s in F1stds]))/4
    F2std = np.sqrt(sum([s**2 for s in F2stds]))/4

    F_lo = (4*(F1+F2-1)+1)/5
    F_lo_std = 4*np.sqrt(F1std**2+F2std**2)/5

    F1, F1std, F2, F2std, F_lo, F_lo_std

    x = list(data.keys())
    y = np.array([data[state]['00']/sum(data[state].values()) for state in x])
    y_shots = F1_shots + F2_shots
    yerr = [np.sqrt(p*(1-p)/y_shots[i]) for i, p in enumerate(y)]

    plt.bar(x, y, yerr=yerr)
    plt.xlabel('Input State', fontsize=12)
    plt.ylabel('State Fidelity', fontsize=12)
    if save == False:
        plt.title('H2-1 logical control-T test, MeasDecode')
    plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    if save == True:
        plt.savefig('control_T_plot_new.pdf', format='pdf')
    plt.show()

    print(f'XZ basis average fidelity: {round(F1,3)} +/- {round(F1std,3)}')
    print(f'ZX basis average fidelity: {round(F2,3)} +/- {round(F2std,3)}')
    print(f'Avg. Fidelity lower bound: {round(F_lo, 3)} +/- {round(F_lo, 3)}')
    
    
def plot_CT_comparison_data(data1, data2, machine=['H2-1', 'H2-1E'], colors=['blue', 'green'], save=False):
    
    F1s = [data1[state]['00']/sum(data1[state].values()) for state in ['+0', '+1', '-0', '-1']]
    F1_shots = [sum(data1[state].values()) for state in ['+0', '+1', '-0', '-1']]
    F1stds = [np.sqrt(p*(1-p)/F1_shots[i]) for i, p in enumerate(F1s)]

    F2s = [data1[state]['00']/sum(data1[state].values()) for state in ['0+', '0-', '1+', '1-']]
    F2_shots = [sum(data1[state].values()) for state in ['0+', '0-', '1+', '1-']]
    F2stds = [np.sqrt(p*(1-p)/F2_shots[i]) for i, p in enumerate(F2s)]
    
    F3s = [data2[state]['00']/sum(data2[state].values()) for state in ['+0', '+1', '-0', '-1']]
    F3_shots = [sum(data2[state].values()) for state in ['+0', '+1', '-0', '-1']]
    F3stds = [np.sqrt(p*(1-p)/F3_shots[i]) for i, p in enumerate(F3s)]

    F4s = [data2[state]['00']/sum(data2[state].values()) for state in ['0+', '0-', '1+', '1-']]
    F4_shots = [sum(data2[state].values()) for state in ['0+', '0-', '1+', '1-']]
    F4stds = [np.sqrt(p*(1-p)/F4_shots[i]) for i, p in enumerate(F4s)]


    F1 = np.mean(F1s)
    F2 = np.mean(F2s)
    
    F3 = np.mean(F3s)
    F4 = np.mean(F4s)

    F1std = np.sqrt(sum([s**2 for s in F1stds]))/4
    F2std = np.sqrt(sum([s**2 for s in F2stds]))/4
    
    F3std = np.sqrt(sum([s**2 for s in F3stds]))/4
    F4std = np.sqrt(sum([s**2 for s in F4stds]))/4

    F_lo1 = (4*(F1+F2-1)+1)/5
    F_lo_std1 = 4*np.sqrt(F1std**2+F2std**2)/5
    
    F_lo2 = (4*(F3+F4-1)+1)/5
    F_lo_std2 = 4*np.sqrt(F3std**2+F4std**2)/5

    #F1, F1std, F2, F2std, F_lo, F_lo_std

    x = list(data1.keys())
    y1_shots = F1_shots + F2_shots
    y2_shots = F3_shots + F4_shots
    y1 = np.array([data1[state]['00']/sum(data1[state].values()) for state in x])
    y2 = np.array([data2[state]['00']/sum(data2[state].values()) for state in x])
    y1err = [np.sqrt(p*(1-p)/y1_shots[i]) for i, p in enumerate(y1)]
    y2err = [np.sqrt(p*(1-p)/y2_shots[i]) for i, p in enumerate(y2)]

    w = 0.4
    plt.bar(np.array(list(range(len(x))))-w/2, y1, width=w, yerr=y1err, label=machine[0], color=colors[0])
    plt.bar(np.array(list(range(len(x))))+w/2, y2, width=w, yerr=y2err, label=machine[1], color=colors[1])
    plt.xlabel('Input State', fontsize=12)
    plt.ylabel('State Fidelity', fontsize=12)
    if save == False:
        plt.title('H2-1 logical control-T test, MeasDecode')
    plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], fontsize=12)
    plt.xticks(ticks=list(range(len(x))), labels=x, fontsize=12)
    plt.legend(loc=1)
    plt.ylim(0,1.1)
    if save == True:
        plt.savefig('control_T_comparison_plot_new.pdf', format='pdf')
    plt.show()
    
    print(machine[0])
    print(f'XZ basis average fidelity: {round(F1,3)} +/- {round(F1std,3)}')
    print(f'ZX basis average fidelity: {round(F2,3)} +/- {round(F2std,3)}')
    print(f'Avg. Fidelity lower bound: {round(F_lo1, 3)} +/- {round(F_lo_std1, 3)}')
    
    print(machine[1])
    print(f'XZ basis average fidelity: {round(F3,3)} +/- {round(F3std,3)}')
    print(f'ZX basis average fidelity: {round(F4,3)} +/- {round(F4std,3)}')
    print(f'Avg. Fidelity lower bound: {round(F_lo2, 3)} +/- {round(F_lo_std2, 3)}')