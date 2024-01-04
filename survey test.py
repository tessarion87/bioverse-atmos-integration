from bioverse.survey import TransitSurvey
from bioverse.generator import Generator
import matplotlib.pyplot as plt
import numpy as np


generator = Generator('transit')
survey = TransitSurvey('default')

t_ref = survey.measurements['has_H2O'].t_ref
T_st_ref = survey.T_st_ref
R_st_ref = survey.R_st_ref


print("Reference star properties: R = {:.1f} R_sun and T_eff = {:.0f} K".format(R_st_ref, T_st_ref))
print("Time required to characterize an Earth-like planet around this star: {:.1f} d".format(t_ref))

sample, detected, data = survey.quickrun(generator, t_total=10*365.25)
t_exp, N_obs = survey.measurements['has_H2O'].compute_exposure_time(data[detected['EEC']])

fig, ax = plt.subplots(ncols=2, figsize=(16,8))

bins = np.logspace(np.log10(0.01), np.log10(np.amax(t_exp)), 30)
ax[0].hist(t_exp, bins=bins)
ax[0].set_xscale('log')
ax[0].set_xlabel('Exposure time (d)')
ax[0].set_ylabel('Number of EECs')
ax[0].axvline(1000/24, linestyle='dashed', lw=5, c='black')

bins = np.logspace(0, 5, 30)
ax[1].hist(N_obs, bins=bins)
ax[1].set_xscale('log')
ax[1].set_xlabel('Number of transit observations')
ax[1].axvline(survey.N_obs_max, linestyle='dashed', lw=5, c='black')

plt.subplots_adjust(wspace=0.3)

plt.show()
plt.clf()

obs = ~np.isnan(data['has_H2O'])
EEC = detected['EEC']

bins = np.linspace(0, 150, 15)
plt.hist(data['d'][obs&EEC], density=True, histtype='step', lw=2, bins=bins, label='Observed')
plt.hist(data['d'][~obs&EEC], density=True, histtype='step', lw=2, bins=bins, label='Not observed')
plt.xlabel('Distance (pc)')
plt.ylabel('Fraction of EECs')
plt.legend(loc='upper left')

plt.show()
