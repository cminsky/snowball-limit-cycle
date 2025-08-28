import numpy as np
import matplotlib.pyplot as plt
import warnings

def savefig(fig,figname,savefmt='pdf',supp_folder=None):
    if supp_folder:
        figname = supp_folder + "/" + figname
    fig.savefig('figs/'+figname,format=savefmt,bbox_inches='tight',dpi=600)

def create_vars(namespace, results_dict):
    for name, array in results_dict.items():
        namespace[name] = array

param_info = {
    # model setup
    'dt': {'unit': 'Myr', 'desc': 'Model time step'},
    't_max': {'unit': 'Myr', 'desc': 'Total model run time'},
    
    # initial conditions
    'T0': {'unit': 'K', 'desc': 'Initial temperature'},
    'a0': {'unit':'','desc':'Initial albedo'},
    
    # ice-albedo feedback parameters
    'Ti': {'unit': 'K', 'desc': 'Ice-covered threshold temperature'},
    'To': {'unit': 'K', 'desc': 'Ice-free threshold temperature'},
    'ai': {'unit': '', 'desc': 'Albedo (ice-covered)'},
    'ao': {'unit': '', 'desc': 'Albedo (ice-free)'},
    
    # background conditions
    'pCO20': {'unit': 'ppm', 'desc': 'Initial atmospheric pCO$_2$'},
    'O20': {'unit': 'Fraction of PAL', 'desc': 'Initial atmospheric O$_2$'},
    
    # volcanism
    'V_C': {'unit': 'Examol/Myr', 'desc': 'Volcanic CO$_2$'},
    'V_red': {'unit': 'Examol/Myr', 'desc': 'Volcanic reduced gases'},
    
    # global weathering
    'W_sil': {'unit': 'Examol/Myr', 'desc': 'Silicate weathering'},
    'W_sea': {'unit': 'Examol/Myr', 'desc': 'Seafloor weathering'},
    'W_org': {'unit':'Examol/Myr', 'desc': 'Organic weathering'},
    'n': {'unit': '', 'desc': 'Silicate weathering feedback strength'},
    
    # biology and burial
    'BorgC': {'unit': 'Examol/Myr', 'desc': 'Organic carbon burial'},
    'forg': {'unit': '', 'desc': 'Fraction of buried C that is organic'},
    'CPsed': {'unit': 'mol/mol', 'desc': 'C:P ratio in sedimentary organic matter'},
    'nb': {'unit': '', 'desc': 'Exponent for org. burial dependence on P'},
    
    # perturbation and LIP parameters
    'C_imb': {'unit': 'Examol/Myr', 'desc': 'Generalized C cycle imbalance'},
    'tau': {'unit': 'Myr', 'desc': 'Timescale of carbon imbalance decay'},
    'W_LIP': {'unit': 'Examol/Myr', 'desc': 'LIP weathering flux'},
    'W_LIP0': {'unit': 'Examol/Myr', 'desc': 'Initial LIP weathering flux'},
    'n_LIP': {'unit': '', 'desc': 'LIP weathering feedback exponent'},
    'PC_LIP': {'unit': '', 'desc': 'P:Ca+Mg molar ratio in LIP rocks'},
    'Wpho_LIP': {'unit': 'Examol/Myr', 'desc': 'Phosphorus weathering from LIP'},
    
    # LIP diagnostics
    'P_eff': {'unit': 'Examol/Myr', 'desc': 'LIP effective weathering'},
    'H_reg': {'unit': 'm', 'desc': 'LIP regolith height'},
    'A_eff': {'unit': '10$^12$ m$^2$', 'desc': 'LIP effective area'},
    'LIP_height':{'unit':'km', 'desc':'LIP total height used up'},
    'LIP_CO2': {'unit': 'Examol', 'desc': 'LIP CO2 drawdown total'},
    'LIP_vol':{'unit': 'm$^3$', 'desc': 'LIP volume used up'},
    
    # phosphorus cycle setup and diagnostics
    'P_conc': {'unit': 'µM', 'desc': 'Seawater phosphate concentration'},
    'W_pho': {'unit': 'Examol/Myr', 'desc': 'Phosphorus weathering'},
    'P0': {'unit': 'Examol', 'desc': 'Initial phosphate reservoir'},
    'P': {'unit': 'Examol', 'desc': 'Phosphate reservoir'},
    'f_Porg': {'unit': '', 'desc': 'Fraction of P buried organically'},
    'Pbur_inorg': {'unit': 'Examol/Myr', 'desc': 'Inorganic phosphorus burial'},
    'Pbur_org': {'unit': 'Examol/Myr', 'desc': 'Org. phosphorus burial'},
    
    # climate diagnostics
    'Teq': {'unit': 'K', 'desc': 'Equilibrium temperature'},
    'a': {'unit': '', 'desc': 'Albedo'},
    'T': {'unit': 'K', 'desc': 'Temperature'},
    'pCO2': {'unit': 'ppm', 'desc': 'pCO$_2$'},
    
    # oxygen cycle diagnostics
    'O2': {'unit': 'Examol', 'desc': 'O2'},
    'pO2': {'unit': 'Fraction of PAL', 'desc': 'pO2 [PAL]'},
    'dOdt': {'unit': 'Examol/Myr', 'desc': 'Oxygen change rate'},
    
    # carbon cycle diagnostics
    'N': {'unit': 'Examol', 'desc': 'Total surficial C'},
    'Nat': {'unit': 'Examol', 'desc': 'Surficial C in atmos'},
    'dNdt': {'unit': 'Examol/Myr', 'desc': 'Change in surficial C change'},
    'dNatdt': {'unit': 'Examol/Myr', 'desc': 'Change in atmospheric C'},
    'F': {'unit': '', 'desc': 'Fraction of surficial C in atmos'},
    'f': {'unit': '', 'desc': 'Fraction of new C stays atmos'},
    'pH': {'unit': '', 'desc': 'pH'},
    
}

from matplotlib import rcParams

rcParams['font.family'] = 'Arial'
rcParams['axes.labelsize'] = 10     # x/y axis labels
rcParams['xtick.labelsize'] = 8     # x-axis tick labels
rcParams['ytick.labelsize'] = 8     # y-axis tick labels
rcParams['legend.fontsize'] = 8     # legend text
rcParams['font.size'] = 8           # default text size
    
# --- plot temperature ---

def annotate_temp(ax,t_max=80,Sturtian=True,Ti=260,To=295):
    ax.set_xlabel('Time (Myr)')
    ax.set_ylabel('Temperature (K)')

    ax.axhspan(To,340, alpha=0.2, color='burlywood', label='Ice-free')
    ax.axhspan(220,Ti, alpha=0.2, color='turquoise', label='Ice-covered')
    ax.set_ylim(230,325)

    xmin, xmax = 0,t_max
    margin = 0.05 * (xmax - xmin)
    ax.set_xlim(xmin - margin, xmax + margin)

    xmax = ax.get_xlim()[1]
    padding = xmax * 0.03
    text_x = xmax - padding

    ax.text(text_x,ax.get_ylim()[0]+5,"ice-covered",c='#10948a',ha='right',va='bottom')
    ax.text(text_x,ax.get_ylim()[1]-5,"ice-free",c='#8b6d32',ha='right',va='top')

    if Sturtian:
        #ax.axvspan(56,59,alpha=0.1,color='k')
        y = 316  # vertical position of the line
        tick_height = 2
        ax.hlines(y, 0, 56, color='k')
        ax.plot([0, 0], [y - tick_height, y + tick_height], c='k')
        ax.plot([56, 56], [y - tick_height, y + tick_height], c='k')
        ax.text(28, y+1, "Sturtian duration", ha='center', va='bottom', fontsize=10)

def plot_temp(results,t_max=80,Sturtian=True,figname=None,savefmt='png'):
    t = results['t']
    T = results['T']
    
    fig, ax = plt.subplots(figsize=(8, 5))
    
    ax.plot(t, T,c='k',alpha=1)
    
    Ti,To = results['Ti'][0],results['To'][0]
    annotate_temp(ax,t_max=t_max,Sturtian=Sturtian,Ti=Ti,To=To)
    
    if figname:
        savefig(fig,figname=figname,savefmt=savefmt)

    plt.show()
    return fig
    
# --- plot oxygen ---

def annotate_pO2(ax,t_max=80):
    ax.set_xlabel('Time (Myr)')
    ax.set_ylabel('pO2 (fraction of PAL)')

    ax.set_ylim(-0.05,0.45)
    
    xmin, xmax = 0,t_max
    margin = 0.05 * (xmax - xmin)
    ax.set_xlim(xmin - margin, xmax + margin)

    xmax = ax.get_xlim()[1]
    padding = xmax * 0.03
    text_x = xmax - padding

    SMIF_threshold = (2.1/1e6)/(0.2) # ppm to % of PAL
    ax.axhspan(ax.get_ylim()[0],SMIF_threshold,color='mediumorchid',alpha=0.3)
    ax.text(text_x,-0.025,'S-MIF expected',color='darkorchid',ha='right',va='center')

    ax.axhspan(0.4,ax.get_ylim()[1],color='g',alpha=0.2)
    ax.text(text_x,0.425,'oxic deep ocean expected',color='darkgreen',ha='right',va='center')
    
def annotate_pO2_yaxis(ax):
    ax.set_ylabel('pO2 (fraction of PAL)')

    ax.set_ylim(-0.05,0.45)

    xmax = ax.get_xlim()[1]
    padding = xmax * 0.03
    text_x = xmax - padding

    SMIF_threshold = (2.1/1e6)/(0.2) # ppm to % of PAL
    ax.axhspan(ax.get_ylim()[0],SMIF_threshold,color='mediumorchid',alpha=0.3)
    ax.text(text_x,-0.025,'S-MIF expected',color='darkorchid',ha='right',va='center')

    ax.axhspan(0.4,ax.get_ylim()[1],color='g',alpha=0.2)
    ax.text(text_x,0.425,'oxic deep ocean expected',color='darkgreen',ha='right',va='center')
    
def plot_O2(results,t_max=80,figname=None,savefmt='png'):    
    t = results['t']
    pO2 = results['pO2']
    
    fig,ax = plt.subplots(figsize=(8, 5))
    
    ax.plot(t, pO2,c='k',alpha=1)
    
    annotate_pO2(ax,t_max=t_max)
    
    if figname:
        savefig(fig,figname=figname,savefmt=savefmt)

    plt.show()
    return fig
    
# --- plot temperature and O2 ---
def plot_main(results,t_max=80,figname=None,savefmt='png',Sturtian=True,vertical=True):
    t = results['t']
    T = results['T']
    pO2 = results['pO2']
    

    if vertical:
        fig, axs = plt.subplots(2,1,figsize=(8,10))
    else:
        fig,axs = plt.subplots(1,2,figsize=(14,5))

    ax = axs[0]
    ax.plot(t, T,c='k',alpha=1)
    Ti,To = results['Ti'][0],results['To'][0]
    annotate_temp(ax,t_max=t_max,Sturtian=Sturtian,Ti=Ti,To=To)

    ax = axs[1]
    ax.plot(t, pO2,c='k',alpha=1)
    annotate_pO2(ax,t_max=t_max)

    if figname:
        savefig(fig,figname=figname,savefmt=savefmt)

    plt.show()
    return fig

# --- analyze time intervals ---

def get_times(t, snowball, return_dict=False):
    snowball_starts = []
    snowball_ends = []
    n = len(snowball)

    for i in range(1, n):
        # start of snowball
        if snowball[i] and not snowball[i-1]:
            snowball_starts.append(t[i])
        # end of snowball
        if not snowball[i] and snowball[i-1]:
            snowball_ends.append(t[i-1])

    # handle case where series ends within a snowball
    if n > 0 and snowball[-1]:
        snowball_ends.append(t[-1])

    # compute snowball durations
    snowball_durations = [end - start for start, end in zip(snowball_starts, snowball_ends)]

    # compute interglacial durations
    interglacial_durations = []
    for end, next_start in zip(snowball_ends, snowball_starts[1:]):
        interglacial_durations.append(next_start - end)
        
    if return_dict:
        return {'snowball_starts':snowball_starts,
               'snowball_ends':snowball_ends,
               'snowball_durations':snowball_durations,
               'interglacial_durations':interglacial_durations}

    return snowball_starts, snowball_ends, snowball_durations, interglacial_durations

def get_end_time(results):
    t = results['t']
    snowball = results['snowball']
    snowball_starts, snowball_ends, snowball_durations, interglacial_durations = get_times(t, snowball)
    end_time = snowball_ends[-1]
    return end_time

def LIP_volume(results,
               Xm=10137, # Ca + Mg mol/m3 basalt
               verbose=True,
               A0=None,
              ):
    W_LIP = np.array(results['W_LIP'])*1e18 # Emol/Myr --> mol/Myr
    dt = results['dt'][0]                    # timestep

    # weathered volume in m3/Myr, integrate over time
    dvol = W_LIP / Xm * dt     # volume increment each step
    vol = np.cumsum(dvol)                # total volume in m^3
    
    if verbose:
        print(f"Total LIP volume used: {vol[-1]/1e15:0.1f} Mkm3")
        
    if A0:
        height = vol/(A0*1e12)
        if verbose:
            print(f"Total LIP height used: {height[-1]/1e3:.2} km")
            P_eff = W_LIP[0]/(A0*1e12)/Xm # m/Myr
            print(f"Initial P_eff required: {P_eff:.0f} m/Myr")
    else:
        height = None
            
    
    return np.array(results['t']), vol, height