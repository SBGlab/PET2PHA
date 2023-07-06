import cobra

#load iJNP4SB_ME
model = cobra.io.load_matlab_model('iJNP4SB_ME.mat')

#Simplification
#STEP 1: Tol plasmid removal
#Tol plasmid genes deletion
plasmidList = ['pWW0_090',
            'pWW0_091',
            'pWW0_092',
            'pWW0_093',
            'pWW0_094',
            'pWW0_095',
            'pWW0_096',
            'pWW0_097',
            'pWW0_099',
            'pWW0_100',
            'pWW0_101',
            'pWW0_102',
            'pWW0_127',
            'pWW0_128',
            'pWW0_129',
            'pWW0_130',
            'pWW0_131']

#PHA Rxns simplification (PHADPC80)
PHAList = [        '3HAACOAT100',
                    '3HAACOAT120',
                    '3HAACOAT121',
                    '3HAACOAT140',
                    '3HAACOAT141',
                    '3HAACOAT60',
                    'ACSPHAC100',
                    'ACSPHAC101',
                    'ACSPHAC120',
                    'ACSPHAC121',
                    'ACSPHAC121d6',
                    'ACSPHAC140',
                    'ACSPHAC141',
                    'ACSPHAC141d5',
                    'ACSPHAC142',
                    'ACSPHAC160',
                    'ACSPHAC40',
                    'ACSPHAC50',
                    'ACSPHAC60',
                    'ACSPHAC70',
                    'ACSPHAC90',
                    'ACSPHACP100',
                    'ACSPHACP40',
                    'ACSPHACP50',
                    'ACSPHACP60',
                    'ACSPHACP70',
                    'ACSPHACP80',
                    'ACSPHACP90',
                    'ACSPHACT40',
                    'ACSPHACT60',
                    'PHADPC100',
                    'PHADPC101',
                    'PHADPC120',
                    'PHADPC121',
                    'PHADPC121d6',
                    'PHADPC140',
                    'PHADPC141',
                    'PHADPC141d5',
                    'PHADPC142',
                    'PHADPC40',
                    'PHADPC50',
                    'PHADPC60',
                    'PHADPC70',
                    'PHADPC90',
                    'PHADPCP100',
                    'PHADPCP40',
                    'PHADPCP50',
                    'PHADPCP60',
                    'PHADPCP70',
                    'PHADPCP80',
                    'PHADPCP90',
                    'PHADPCT40',
                    'PHADPCT60',
                    'PHAP2C100',
                    'PHAP2C101',
                    'PHAP2C120',
                    'PHAP2C121',
                    'PHAP2C121d6',
                    'PHAP2C140',
                    'PHAP2C141',
                    'PHAP2C141d5',
                    'PHAP2C142',
                    'PHAP2C60',
                    'PHAPC40',
                    'PHAPC50',
                    'PHAPC70',
                    'PHAPC90',
                    'PHAPCP100',
                    'PHAPCP40',
                    'PHAPCP50',
                    'PHAPCP60',
                    'PHAPCP70',
                    'PHAPCP80',
                    'PHAPCP90',
                    'PHAPCT40',
                    'PHAPCT60',
                    'RECOAH10',
                    'RECOAH11',
                    'RECOAH12',
                    'RECOAH13',
                    'RECOAH14',
                    'RECOAH15',
                    'RECOAH16',
                    'RECOAH17',
                    'RECOAH19',
                    'RECOAH2',
                    'RECOAH20',
                    'RECOAH21',
                    'RECOAH22',
                    'RECOAH23',
                    'RECOAH24',
                    'RECOAH25',
                    'RECOAH26',
                    'RECOAH4',
                    'RECOAH5',
                    'RECOAH6',
                    'RECOAH7',
                    'RECOAH8',
                    'RECOAH9',
                    'RHACOAE100',
                    'RHACOAE120',
                    'RHACOAE121',
                    'RHACOAE121d6',
                    'RHACOAE140',
                    'RHACOAE141',
                    'RHACOAE141d5',
                    'RHACOAE142',
                    'RHACOAE160',
                    'RHACOAE40',
                    'RHACOAE50',
                    'RHACOAE60',
                    'RHACOAE70',
                    'RHACOAE90',
                    'RHACOAEP100',
                    'RHACOAEP101',
                    'RHACOAEP40',
                    'RHACOAEP50',
                    'RHACOAEP60',
                    'RHACOAEP70',
                    'RHACOAEP80',
                    'RHACOAEP90',
                    'RHACOAET40',
                    'RHACOAET60',
                    'RHACOAR100',
                    'RHACOAR101',
                    'RHACOAR120',
                    'RHACOAR121d5',
                    'RHACOAR121d6',
                    'RHACOAR140',
                    'RHACOAR141d5',
                    'RHACOAR141d7',
                    'RHACOAR142',
                    'RHACOAR60',
                    'RHACOAR70',
                    'RHACOAR90',
                    'RHACOARP100',
                    'RHACOARP40',
                    'RHACOARP50',
                    'RHACOARP60',
                    'RHACOARP70',
                    'RHACOARP80',
                    'RHACOARP90',
                    'RHACOART60',
                    'RHA100tpp',
                    'RHA101tpp',
                    'RHA120tpp',
                    'RHA121d6tpp',
                    'RHA121tpp',
                    'RHA140tpp',
                    'RHA141d5tpp',
                    'RHA141tpp',
                    'RHA142tpp',
                    'RHA160tpp',
                    'RHA40tpp',
                    'RHA50tpp',
                    'RHA60tpp',
                    'RHA70tpp',
                    'RHA90tpp',
                    'RHAP100tpp',
                    'RHAP40tpp',
                    'RHAP50tpp',
                    'RHAP60tpp',
                    'RHAP70tpp',
                    'RHAP80tpp',
                    'RHAP90tpp',
                    'RHAT40tpp',
                    'RHAT60tpp',
                    'RHA100tex',
                    'RHA101tex',
                    'RHA120tex',
                    'RHA121d6tex',
                    'RHA121tex',
                    'RHA140tex',
                    'RHA141d5tex',
                    'RHA141tex',
                    'RHA142tex',
                    'RHA160tex',
                    'RHA40tex',
                    'RHA50tex',
                    'RHA60tex',
                    'RHA70tex',
                    'RHA90tex',
                    'RHAP100tex',
                    'RHAP40tex',
                    'RHAP50tex',
                    'RHAP60tex',
                    'RHAP70tex',
                    'RHAP80tex',
                    'RHAP90tex',
                    'RHAT40tex',
                    'RHAT60tex',
                    'DM_C100aPHA',
                    'DM_C100pPHA',
                    'DM_C101PHA',
                    'DM_C120aPHA',
                    'DM_C121aPHA',
                    'DM_C121d6PHA',
                    'DM_C140aPHA',
                    'DM_C141aPHA',
                    'DM_C141d5PHA',
                    'DM_C142PHA',
                    'DM_C40aPHA',
                    'DM_C40atPHA',
                    'DM_C40pPHA',
                    'DM_C50aPHA',
                    'DM_C50pPHA',
                    'DM_C60aPHA',
                    'DM_C60atPHA',
                    'DM_C60pPHA',
                    'DM_C70aPHA',
                    'DM_C70pPHA',
                    'DM_C80pPHA',
                    'DM_C90aPHA',
                    'DM_C90pPHA']


#STEP3: PDH, AKGDH Rxns deletion
DHList = ['PDH','AKGDH']


#STEP4: Alginate Rxns simplification (23)
AlginateList = ['ALGAC5',
                'ALGE3',
                'ALGE5',
                'ALGSex3',
                'ALGSex5',
                'ALGSex6',
                'EX_algac_M_(e)',
                'EX_glc(e)',
                'EX_prealg_MG_23_(e)',
                'EX_prealginate_G_(e)',
                'ALGE1',
                'ALGE2',
                'ALGE4',
                'ALGE6',
                'ALGE7',
                'ALGE9',
                'ALGSex1',
                'ALGSex10',
                'ALGSex2',
                'ALGSex4',
                'ALGSex7',
                'ALGSex8',
                'EX_algac_MG_14_(e)',
                'EX_algac_MG_32_(e)',
                'EX_algac_MG_41_(e)',
                'EX_prealg_MG_14_(e)',
                'EX_prealg_MG_32_(e)',
                'EX_prealg_MG_41_(e)']

#STEP5: CHANGE BOUNDS TO KEY REACTIONS
knoledge_based_kos = [  'GLCtex', 'TOLt5', 'TOLt6', 'TOLtex']
pre_gdls_kos = ['ALDD2x','CYSS','NACODA','PHADPC80','RHACOAE80']
# are part of the gdls strategy, among them CYSS and NACODA are essential as the removal of one of them
#decouples the production

#Remove the rest of biomass reactions (add to kt_simplification.py)
biomass_rxns = ['BiomassKT2440_Core2','BiomassKT2440_WT3']

def main():
    #DELETE ALL LISTED REACTIONS AND GENES:
    r_kos = AlginateList + DHList + PHAList + knoledge_based_kos + pre_gdls_kos + biomass_rxns
    g_kos =plasmidList
    for ko in r_kos:
        model.reactions.get_by_id(ko).bounds = (0,0)
        
    for ko in g_kos:
        model.genes.get_by_id(ko).knock_out()
        
    #BLOCK GLUCOSE UPTAKE
    glucose_uptake = 'EX_glc(e)'
    model.reactions.get_by_id(glucose_uptake).bounds = (0, 1000)
    
    #Save simplified model:
    cobra.io.write_sbml_model(model, 'iJNP4SB_ME_simplified.xml')


if __name__ == '__main__':
    main()

#SBML model generated with this procedure has an uptake of glucose, ask if it is that way or should be changed to 0
