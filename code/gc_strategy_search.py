#!/usr/bin/env python
# coding: utf-8

# # Strain-desing Workflow

## IMPORTS
#Cobra dependencies
import cobra
from cobra import Reaction, Metabolite
from cobra.io import read_sbml_model, save_matlab_model
from cobra.flux_analysis import production_envelope, flux_variability_analysis
from cobra.flux_analysis.variability import find_essential_genes
from cobra.sampling import sample
from cobra.sampling import OptGPSampler, ACHRSampler
import gurobipy
#Cameo dependencies
import cameo
from cameo.strain_design.deterministic.linear_programming import OptKnock
from cameo import phenotypic_phase_plane
from cameo.visualization.plotting.with_plotly import PlotlyPlotter
from cameo.strain_design.deterministic.flux_variability_based import FSEOF
from cameo.flux_analysis.simulation import lmoma, pfba
#Dependencies for Minimal Cut Sets (MCS) analysis
import straindesign as sd
#Data processing dependencies
import ast
import re
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from itertools import combinations
#Metabolic design helper functions
from pyfastcore import set_medium
#Our own helper functions
import subprocess #to call matlab gcFront execution script
from utils.designFunctions import *  #our own design functions
from utils.importExcelModel import * #to generate GEMs from excel spreadsheets
#widgets
import ipywidgets as widgets
from ipywidgets import Layout, HBox, VBox, HTML
#warnings
import warnings
#For show/hide code
from IPython.display import HTML as hd
#For OS interaction
import os
from os import path
import time

###################################################################################
# 1. SET PARAMETERS

#Set the parameters for exploration of all combinations:                          
max_knock_out_range =  range(3, 6)
max_cl_range = range(9, 13)
replicates = 3

#Project-specific parameters
dir_path = '.' # for actual directory
project_name = 'putida_pet_to_pha_ko_results'
fig_dir_path = '../figures' 
model_path = '../models/iJNP4SB_ME_simplified.xml'
framework_name = 'cameo'

#Bioprocess-specific parameters
carbon_source = 'EX_pet(e)'
uptake = -3.78
target_biomass = 'BiomassKT2440_ME'
target_reaction = 'DM_C80aPHA'
minimum_growth_fraction = 0.10 
min_production_flux = 0.05

#Strategy sorting parameters
sorting_param = 'Normalised_max_Flux'   #TO CHOOSE BETWEEN: 'Normalised_max_Flux',               
secondary_param = 'Presence'            #'Normalised_max_biomass' or 'Presence'

set_up_params = {'dir_path' : dir_path,
                 'project_name' : project_name,
                 'framework' : framework_name,
                 'min_growth' : minimum_growth_fraction,
                 'max_cl_range' : max_cl_range,
                 'max_knock_out_range' : max_knock_out_range,
                 'replicates' : replicates}

###################################################################################

# COMPUTATIONAL WORKFLOW


def main():
# 2. Configure the model according to bioprocess conditions

# 2.1 Model & media specification

    # Specify the model and the required parameters
    model=cobra.io.read_sbml_model(model_path)
    target_metabolites = [m.id for m in model.reactions.get_by_id(target_reaction).products] # for the example ['btd_RR_e']
    #m9 media as default media
    media_definition = {'EX_ca2(e)' : -1000,
                        'EX_cl(e)' : -1000,
                        'EX_co2(e)' : -1000,
                        'EX_cobalt2(e)' : -1000,
                        'EX_cu2(e)' : -1000,
                        'EX_fe2(e)' : -1000,
                        'EX_fe3(e)' : -1000,
                        'EX_h(e)' : -1000,
                        'EX_h2o(e)' : -1000,
                        'EX_k(e)' : -1000,
                        'EX_mg2(e)' : -1000,
                        'EX_mn2(e)' : -1000,
                        'EX_mobd(e)' : -1000,
                        'EX_na1(e)' : -1000,
                        'EX_tungs(e)' : -1000,
                        'EX_zn2(e)' : -1000,
                        'EX_ni2(e)' : -1000,
                        'EX_sel(e)' : -1000,
                        'EX_slnt(e)' : -1000,
                        'EX_so4(e)' : -1000,
                        'EX_nh4(e)' : -1000,
                        'EX_pi(e)' : -1000,
                        'EX_cbl1(e)' : -.01,
                        'EX_o2(e)' : -20,
                        carbon_source : uptake
                       }

    model = set_medium(model, media_definition)

# 2.2 Model reduction

    # Preprocess the model with user data
    configured_model = model.copy()
    configured_model.objective = target_biomass
    # Specify media if not previously define in a preparation script
    configured_model = purge_non_objective_biomass(configured_model, target_biomass, n_of_biomass_reactions=3)
    blocked_reactions = set([r.id for r in configured_model.reactions])-set(get_rxn_with_fva_flux(configured_model))
    #remove blocked reactions
    print('Removing a total of %d blocked reactions...' % len(blocked_reactions))
    configured_model.remove_reactions(blocked_reactions)
    # Save model to .mat format if it were neccessary to run gcFront after
    save_matlab_model(configured_model, model_path.replace("iJNP4SB_ME_simplified.xml", "configured_model.mat") )
    wt_growth = configured_model.optimize().fluxes[target_biomass]
    biomas_limit = wt_growth*minimum_growth_fraction


# 3.Generate a candidate reaction list through OptKnock

    #FIRST: we execute several optknock runs over the parameter space for all the specified replicates:
    print('EXECUTING OPTKNOCK SEARCH OF STRATEGIES...')
    analysis_type = ['gene_candidates','strategies_eval']
    generate_strategies(set_up_params, configured_model, target_biomass, target_reaction, carbon_source, blocked_reactions)   
    #SECOND: analysis of OptKnock results
    print('ANALYSIS OF OPTKNOCK RESULTS')
    ok_raw_df, ok_results_df = analyse_results(set_up_params, configured_model, target_reaction, sorting_param, analysis_type = analysis_type)

    #ranking figure of top reactions in OptKnock result:
    print('Ranking top reactions in Optknock designs...')
    #delete all reactions with a sorting_param lower than 10% of best strategy
    param_max = max(ok_results_df[sorting_param].tolist())
    ok_results_df = ok_results_df.loc[ok_results_df[sorting_param]>=0.1*param_max]
    ranking_fig = px.bar(ok_results_df, x="Reaction", y=sorting_param, color=secondary_param, color_continuous_scale='Bluered_r')
    
    #strategies performance over the parameter space:
    deletion_list = ok_raw_df.query("Flux>="+str(min_production_flux))["N_of_deletions"].unique().tolist()
    deletion_list.sort()
    carbon_list = ok_raw_df.query("Flux>="+str(min_production_flux))["C_limit"].unique().tolist()
    carbon_list.sort()
    performance_fig = px.scatter(ok_raw_df.query("Flux>="+str(min_production_flux)),
                                 x="Flux", y="Biomass", color="Biomass",
                                 facet_col="N_of_deletions", facet_row="C_limit",hover_name="Strategy",
                                 category_orders = {"N_of_deletions":deletion_list,
                                                    "C_limit":carbon_list})

    performance_fig.for_each_annotation(lambda a: a.update(text=a.text.replace("N_of_deletions", "dels")))

    print('UserWarning : Sometimes Optknock strategies are not reproducible when applied over the model, please check if they work')

    #save figures:
    ranking_fig.write_image('/'.join([fig_dir_path, "candidate_reactions_reduced_parameter_space.svg"]))
    performance_fig.write_image('/'.join([fig_dir_path, "optknock_strategies_evaluation_reduced_parameter_space.svg"]))

# 4. Execute the suitable strategy design algortihm

    #Now the workflow branch for searching gc strategies should start
    #Execution of gcFront
    print('EXECUTING gcFront in MATLAB...')
    #Generate a file containing basic information of the set up
    set_up_dict = {'experiment' : [project_name],
                   'target' : [target_reaction],
                   'framework' : [framework_name],
                   'biomass_limit' : [biomas_limit],
                   'minimum_flux' : [min_production_flux]
                  }

    set_up_df = pd.DataFrame.from_dict(set_up_dict)
    set_up_df.to_csv('production_target_tf.csv',index=False)
    #Call bash to launch the gcFront script in the background
    subprocess.call("bash utils/launch_gcFront.sh", shell=True)
    print("gcFront was executed, the process could take hours, meanwhile you can explore OptKnock strategies")
    print("We reccomend you to wait for results and analyse strategies in 'Strategies_Analysis.ipynb' notebook")

# END

if __name__ == '__main__':
    main()
    
