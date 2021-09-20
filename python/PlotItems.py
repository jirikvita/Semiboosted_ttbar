#!/usr/bin/python

LastCut = 'AnySel'
FirstCut = 'NoCuts'
TagDir = 'TagHistos'
migradir = 'migrations/'

# suffix for more dense version of plots, 4.8.2020
suff = '_denser'
nosuff = ''

names1d_triple =[     
#[LastCut + '/Detector_partonDijetMass', LastCut + '/ParticleDijetMass'] #, LastCut + '/DetectorDijetMass'],
#    [LastCut + '/Detector_partonDijetMass', LastCut + '/ParticleDijetMass'] 
    [LastCut + '/Detector_partonDijetMass', LastCut + '/DetectorDijetMass'],
]


# TODO:
# make a generation script;)

names1d_Glob = [

    # DiLead stuff, 2.7.2021
    #['AnySel/LeadDiParticleJetdR', 'AnySel/LeadDiDetectorJetdR'],
    #['AnySel/LeadDiParticleJetMass', 'AnySel/LeadDiDetectorJetMass'],
    #['AnySel/LeadDiParticlebJetdR', 'AnySel/LeadDiDetectorbJetdR'],
    #['AnySel/LeadDiParticlebJetMass', 'AnySel/LeadDiDetectorbJetMass'],
    
    # globals:
    ['AnySel/ParticleMet', 'AnySel/DetectorMet'],
    ['AnySel/ParticleHTj', 'AnySel/DetectorHTj'],
    ['AnySel/ParticleHTJ', 'AnySel/DetectorHTJ'],
    ['AnySel/ParticleHTjPlusMet', 'AnySel/DetectorHTjPlusMet'],
    ['AnySel/ParticleHTJPlusMet', 'AnySel/DetectorHTJPlusMet'],
    ['AnySel/ParticleSumMj', 'AnySel/DetectorSumMj'],
    ['AnySel/ParticleSumMJ', 'AnySel/DetectorSumMJ'],
    ['AnySel/ParticleMjVis', 'AnySel/DetectorMjVis'],
    ['AnySel/ParticleMJVis', 'AnySel/DetectorMJVis'],
    ['AnySel/ParticlejAplanarity', 'AnySel/DetectorjAplanarity'],
    ['AnySel/ParticlejSphericity', 'AnySel/DetectorjSphericity'],
    ['AnySel/ParticleJAplanarity', 'AnySel/DetectorJAplanarity'],
    ['AnySel/ParticleJSphericity', 'AnySel/DetectorJSphericity'],

    ['2B0S/ParticleMet', '2B0S/DetectorMet'],
    ['2B0S/ParticleHTj', '2B0S/DetectorHTj'],
    ['2B0S/ParticleHTJ', '2B0S/DetectorHTJ'],
    ['2B0S/ParticleHTjPlusMet', '2B0S/DetectorHTjPlusMet'],
    ['2B0S/ParticleHTJPlusMet', '2B0S/DetectorHTJPlusMet'],
    ['2B0S/ParticleSumMj', '2B0S/DetectorSumMj'],
    ['2B0S/ParticleSumMJ', '2B0S/DetectorSumMJ'],
    ['2B0S/ParticleMjVis', '2B0S/DetectorMjVis'],
    ['2B0S/ParticleMJVis', '2B0S/DetectorMJVis'],    
    ['2B0S/ParticlejAplanarity', '2B0S/DetectorjAplanarity'],
    ['2B0S/ParticlejSphericity', '2B0S/DetectorjSphericity'],
    ['2B0S/ParticleJAplanarity', '2B0S/DetectorJAplanarity'],
    ['2B0S/ParticleJSphericity', '2B0S/DetectorJSphericity'],

    ['1B1S/ParticleMet', '1B1S/DetectorMet'],
    ['1B1S/ParticleHTj', '1B1S/DetectorHTj'],
    ['1B1S/ParticleHTJ', '1B1S/DetectorHTJ'],
    ['1B1S/ParticleHTjPlusMet', '1B1S/DetectorHTjPlusMet'],
    ['1B1S/ParticleHTJPlusMet', '1B1S/DetectorHTJPlusMet'],
    ['1B1S/ParticleSumMj', '1B1S/DetectorSumMj'],
    ['1B1S/ParticleSumMJ', '1B1S/DetectorSumMJ'],
    ['1B1S/ParticleMjVis', '1B1S/DetectorMjVis'],
    ['1B1S/ParticleMJVis', '1B1S/DetectorMJVis'],    
    ['1B1S/ParticlejAplanarity', '1B1S/DetectorjAplanarity'],
    ['1B1S/ParticlejSphericity', '1B1S/DetectorjSphericity'],
    ['1B1S/ParticleJAplanarity', '1B1S/DetectorJAplanarity'],
    ['1B1S/ParticleJSphericity', '1B1S/DetectorJSphericity'],

    ['0B2S/ParticleMet', '0B2S/DetectorMet'],
    ['0B2S/ParticleHTj', '0B2S/DetectorHTj'],
    ['0B2S/ParticleHTJ', '0B2S/DetectorHTJ'],
    ['0B2S/ParticleHTjPlusMet', '0B2S/DetectorHTjPlusMet'],
    ['0B2S/ParticleHTJPlusMet', '0B2S/DetectorHTJPlusMet'],
    ['0B2S/ParticleSumMj', '0B2S/DetectorSumMj'],
    ['0B2S/ParticleSumMJ', '0B2S/DetectorSumMJ'],
    ['0B2S/ParticleMjVis', '0B2S/DetectorMjVis'],
    ['0B2S/ParticleMJVis', '0B2S/DetectorMJVis'],    
    ['0B2S/ParticlejAplanarity', '0B2S/DetectorjAplanarity'],
    ['0B2S/ParticlejSphericity', '0B2S/DetectorjSphericity'],
    ['0B2S/ParticleJAplanarity', '0B2S/DetectorJAplanarity'],
    ['0B2S/ParticleJSphericity', '0B2S/DetectorJSphericity'],
]

names1d_Rel = [

       
    # relative variables: 28.5.2020
    ['2B0S/ParticleDiTopPtRel',   '2B0S/DetectorDiTopPtRel'],
    ['2B0S/ParticleDiTopPtGeo',   '2B0S/DetectorDiTopPtGeo'],
    ['2B0S/ParticleDiTopMassGeo',   '2B0S/DetectorDiTopMassGeo'],
    ['2B0S/ParticleTopPtRel',   '2B0S/DetectorTopPtRel'],
    ['2B0S/ParticleTop1PtRel',   '2B0S/DetectorTop1PtRel'],
    ['2B0S/ParticleTop2PtRel',   '2B0S/DetectorTop2PtRel'],
    ['2B0S/ParticleDiTopPoutRel',   '2B0S/DetectorDiTopPoutRel'],
    ['2B0S/ParticleDiTopPoutGeo',   '2B0S/DetectorDiTopPoutGeo'],

    ['1B1S/ParticleDiTopPtRel',   '1B1S/DetectorDiTopPtRel'],
    ['1B1S/ParticleDiTopPtGeo',   '1B1S/DetectorDiTopPtGeo'],
    ['1B1S/ParticleDiTopMassGeo',   '1B1S/DetectorDiTopMassGeo'],
    ['1B1S/ParticleTopPtRel',   '1B1S/DetectorTopPtRel'],
    ['1B1S/ParticleTop1PtRel',   '1B1S/DetectorTop1PtRel'],
    ['1B1S/ParticleTop2PtRel',   '1B1S/DetectorTop2PtRel'],
    ['1B1S/ParticleDiTopPoutRel',   '1B1S/DetectorDiTopPoutRel'],
    ['1B1S/ParticleDiTopPoutGeo',   '1B1S/DetectorDiTopPoutGeo'],

    ['0B2S/ParticleDiTopPtRel',   '0B2S/DetectorDiTopPtRel'],
    ['0B2S/ParticleDiTopPtGeo',   '0B2S/DetectorDiTopPtGeo'],
    ['0B2S/ParticleDiTopMassGeo',   '0B2S/DetectorDiTopMassGeo'],
    ['0B2S/ParticleTopPtRel',   '0B2S/DetectorTopPtRel'],
    ['0B2S/ParticleTop1PtRel',   '0B2S/DetectorTop1PtRel'],
    ['0B2S/ParticleTop2PtRel',   '0B2S/DetectorTop2PtRel'],
    ['0B2S/ParticleDiTopPoutRel',   '0B2S/DetectorDiTopPoutRel'],
    ['0B2S/ParticleDiTopPoutGeo',   '0B2S/DetectorDiTopPoutGeo']
]


names1d_SBreco = [


    # masses:
    
    ['2B0S/ParticleTop1Mass',   '2B0S/DetectorTop1Mass'],
    ['2B0S/ParticleTop2Mass',   '2B0S/DetectorTop2Mass'], 
    ['2B0S/ParticleTopMass',    '2B0S/DetectorTopMass'],

    ['1B1S/ParticleTop1Mass',   '1B1S/DetectorTop1Mass'],
    ['1B1S/ParticleTop2Mass',   '1B1S/DetectorTop2Mass'], 
    ['1B1S/ParticleTopMass',    '1B1S/DetectorTopMass'],

    ['0B2S/ParticleTop1Mass',   '0B2S/DetectorTop1Mass'],
    ['0B2S/ParticleTop2Mass',   '0B2S/DetectorTop2Mass'], 
    ['0B2S/ParticleTopMass',    '0B2S/DetectorTopMass'],
    
    ['1B1S/ParticleW1Mass',   '1B1S/DetectorW1Mass'],
    ['1B1S/ParticleW2Mass',   '1B1S/DetectorW2Mass'], 
    ['1B1S/ParticleWMass',    '1B1S/DetectorWMass'],
    
    ['0B2S/ParticleW1Mass',   '0B2S/DetectorW1Mass'],
    ['0B2S/ParticleW2Mass',   '0B2S/DetectorW2Mass'], 
    ['0B2S/ParticleWMass',    '0B2S/DetectorWMass'],


    # pTs:
    
    ['2B0S/ParticleTop1Pt' + suff,   '2B0S/DetectorTop1Pt' + suff],
    ['2B0S/ParticleTop2Pt' + suff,   '2B0S/DetectorTop2Pt' + suff], 
    ['2B0S/ParticleTopPt' + suff,    '2B0S/DetectorTopPt' + suff],

    ['1B1S/ParticleTop1Pt' + suff,   '1B1S/DetectorTop1Pt' + suff],
    ['1B1S/ParticleTop2Pt' + suff,   '1B1S/DetectorTop2Pt' + suff], 
    ['1B1S/ParticleTopPt' + suff,    '1B1S/DetectorTopPt' + suff],

    ['0B2S/ParticleTop1Pt' + suff,   '0B2S/DetectorTop1Pt' + suff],
    ['0B2S/ParticleTop2Pt' + suff,   '0B2S/DetectorTop2Pt' + suff], 
    ['0B2S/ParticleTopPt' + suff,    '0B2S/DetectorTopPt' + suff],
    
    ['1B1S/ParticleW1Pt',   '1B1S/DetectorW1Pt'],
    ['1B1S/ParticleW2Pt',   '1B1S/DetectorW2Pt'], 
    ['1B1S/ParticleWPt',    '1B1S/DetectorWPt'],
    
    ['0B2S/ParticleW1Pt',   '0B2S/DetectorW1Pt'],
    ['0B2S/ParticleW2Pt',   '0B2S/DetectorW2Pt'], 
    ['0B2S/ParticleWPt',    '0B2S/DetectorWPt'],

    ['2B0S/ParticleDiTopRttbar' + suff,    '2B0S/DetectorDiTopRttbar' + suff],
    #['2B1S/ParticleDiTopRttbar' + suff,    '2B1S/DetectorDiTopRttbar' + suff],
    ['1B1S/ParticleDiTopRttbar' + suff,    '1B1S/DetectorDiTopRttbar' + suff],
    ['0B2S/ParticleDiTopRttbar' + suff,    '0B2S/DetectorDiTopRttbar' + suff],
    
    ['2B0S/ParticleDiTopMass' + nosuff,    '2B0S/DetectorDiTopMass' + nosuff],
    #['2B1S/ParticleDiTopMass' + nosuff,    '2B1S/DetectorDiTopMass' + nosuff],
    ['1B1S/ParticleDiTopMass' + nosuff,    '1B1S/DetectorDiTopMass' + nosuff],
    ['0B2S/ParticleDiTopMass' + nosuff,    '0B2S/DetectorDiTopMass' + nosuff],

    ['2B0S/ParticleDiTopMass' + suff,    '2B0S/DetectorDiTopMass' + suff],
    #['2B1S/ParticleDiTopMass' + suff,    '2B1S/DetectorDiTopMass' + suff],
    ['1B1S/ParticleDiTopMass' + suff,    '1B1S/DetectorDiTopMass' + suff],
    ['0B2S/ParticleDiTopMass' + suff,    '0B2S/DetectorDiTopMass' + suff],

    
    ['2B0S/ParticleDiTopPt' + suff,    '2B0S/DetectorDiTopPt' + suff],
    #['2B1S/ParticleDiTopPt' + suff,    '2B1S/DetectorDiTopPt' + suff],
    ['1B1S/ParticleDiTopPt' + suff,    '1B1S/DetectorDiTopPt' + suff],
    ['0B2S/ParticleDiTopPt' + suff,    '0B2S/DetectorDiTopPt' + suff],



    ['2B0S/ParticleDiTopPout' + suff,    '2B0S/DetectorDiTopPout' + suff],
    #['2B1S/ParticleDiTopPout' + suff,    '2B1S/DetectorDiTopPout' + suff],
    ['1B1S/ParticleDiTopPout' + suff,    '1B1S/DetectorDiTopPout' + suff],
    ['0B2S/ParticleDiTopPout' + suff,    '0B2S/DetectorDiTopPout' + suff],

    ['2B0S/ParticleDiTopChittbar',    '2B0S/DetectorDiTopChittbar'],
    #['2B1S/ParticleDiTopChittbar',    '2B1S/DetectorDiTopChittbar'],
    ['1B1S/ParticleDiTopChittbar',    '1B1S/DetectorDiTopChittbar'],
    ['0B2S/ParticleDiTopChittbar',    '0B2S/DetectorDiTopChittbar'],

    ['2B0S/ParticleDiTopCosThetaStar' + suff,    '2B0S/DetectorDiTopCosThetaStar' + suff],
    #['2B1S/ParticleDiTopCosThetaStar' + suff,    '2B1S/DetectorDiTopCosThetaStar' + suff],
    ['1B1S/ParticleDiTopCosThetaStar' + suff,    '1B1S/DetectorDiTopCosThetaStar' + suff],
    ['0B2S/ParticleDiTopCosThetaStar' + suff,    '0B2S/DetectorDiTopCosThetaStar' + suff],

    ['2B0S/ParticleDiTopDeltaPhi',    '2B0S/DetectorDiTopDeltaPhi'],
    #['2B1S/ParticleDiTopDeltaPhi',    '2B1S/DetectorDiTopDeltaPhi'],
    ['1B1S/ParticleDiTopDeltaPhi',    '1B1S/DetectorDiTopDeltaPhi'],
    ['0B2S/ParticleDiTopDeltaPhi',    '0B2S/DetectorDiTopDeltaPhi'],
    
    ['2B0S/ParticleDiTopDelta',    '2B0S/DetectorDiTopDelta'],
    #['2B1S/ParticleDiTopDelta',    '2B1S/DetectorDiTopDelta'],
    ['1B1S/ParticleDiTopDelta',    '1B1S/DetectorDiTopDelta'],
    ['0B2S/ParticleDiTopDelta',    '0B2S/DetectorDiTopDelta'],

    ['2B0S/ParticleDiTopYboost',    '2B0S/DetectorDiTopYboost'],
    #['2B1S/ParticleDiTopYboost',    '2B1S/DetectorDiTopYboost'],
    ['1B1S/ParticleDiTopYboost',    '1B1S/DetectorDiTopYboost'],
    ['0B2S/ParticleDiTopYboost',    '0B2S/DetectorDiTopYboost'],

    
    # any selection:

    ['AnySel/ParticleDiTopMass' + suff,    'AnySel/DetectorDiTopMass' + suff],
    ['AnySel/ParticleDiTopPt' + suff,    'AnySel/DetectorDiTopPt' + suff],

    ['AnySel/ParticleTop1Mass',   'AnySel/DetectorTop1Mass'],
    ['AnySel/ParticleTop2Mass',   'AnySel/DetectorTop2Mass'], 
    ['AnySel/ParticleTopMass',    'AnySel/DetectorTopMass'],

    ['AnySel/ParticleTop1Pt',   'AnySel/DetectorTop1Pt'],
    ['AnySel/ParticleTop2Pt',   'AnySel/DetectorTop2Pt'], 
    ['AnySel/ParticleTopPt',    'AnySel/DetectorTopPt'],

    ['AnySel/ParticleW1Pt',   'AnySel/DetectorW1Pt'],
    ['AnySel/ParticleW2Pt',   'AnySel/DetectorW2Pt'], 
    ['AnySel/ParticleWPt',    'AnySel/DetectorWPt'],

    


    
]

names1d_gen = [       
        ]


names2d_tagging = [
    
    FirstCut + '/DetectorLJetTau32VsMass',
    FirstCut + '/ParticleLJetTau32VsMass',
    FirstCut + '/DetectorLJetTau21VsMass',
    FirstCut + '/ParticleLJetTau21VsMass',
    FirstCut + '/DetectorLJetTau21VsTau32',
    FirstCut + '/ParticleLJetTau21VsTau32',
    FirstCut + '/DetectorLJetTau32VsPt',
    FirstCut + '/ParticleLJetTau32VsPt',
    FirstCut + '/DetectorLJetTau21VsPt',
    FirstCut + '/ParticleLJetTau21VsPt',
    FirstCut + '/DetectorLJetPtVsMass',
    FirstCut + '/ParticleLJetPtVsMass',
   
    TagDir + '/ttruthMatchedLjetsTau32VsMassParticle',
    TagDir + '/ttruthMatchedLjetsTau32VsMassDetector',
    TagDir + '/WtruthMatchedLjetsTau32VsMassParticle',
    TagDir + '/WtruthMatchedLjetsTau32VsMassDetector',
    TagDir + '/WTagLjetsTau32VsMassDetector',
    TagDir + '/TopTagLjetsTau32VsMassDetector',

    TagDir + '/ttruthMatchedLjetsTau21VsMassParticle',
    TagDir + '/ttruthMatchedLjetsTau21VsMassDetector',
    TagDir + '/WtruthMatchedLjetsTau21VsMassParticle',
    TagDir + '/WtruthMatchedLjetsTau21VsMassDetector',
    TagDir + '/WTagLjetsTau21VsMassDetector',
    TagDir + '/TopTagLjetsTau21VsMassDetector',

    TagDir + '/ttruthMatchedLjetsTau32VsPtParticle',
    TagDir + '/ttruthMatchedLjetsTau32VsPtDetector',
    TagDir + '/WtruthMatchedLjetsTau32VsPtParticle',
    TagDir + '/WtruthMatchedLjetsTau32VsPtDetector',
    TagDir + '/TopTagLjetsTau32VsPtDetector',
    TagDir + '/WTagLjetsTau32VsPtDetector',

    TagDir + '/ttruthMatchedLjetsTau21VsPtParticle',
    TagDir + '/ttruthMatchedLjetsTau21VsPtDetector',
    TagDir + '/WtruthMatchedLjetsTau21VsPtParticle',
    TagDir + '/WtruthMatchedLjetsTau21VsPtDetector',
    TagDir + '/WTagLjetsTau21VsPtDetector',
    TagDir + '/TopTagLjetsTau21VsPtDetector',

    TagDir + '/ttruthMatchedLjetsPtVsMassParticle',
    TagDir + '/ttruthMatchedLjetsPtVsMassDetector',
    TagDir + '/WtruthMatchedLjetsPtVsMassParticle',
    TagDir + '/WtruthMatchedLjetsPtVsMassDetector',
    TagDir + '/WTagLjetsPtVsMassDetector',
    TagDir + '/TopTagLjetsPtVsMassDetector',

    TagDir +'/LjetPtVsMinDRLjetWPartonZoomDetector',
    TagDir +'/LjetPtVsMinDRLjetWPartonZoomParticle',
    TagDir +'/LjetPtVsMinDRLjetWPartonDetector',
    TagDir +'/LjetPtVsMinDRLjetWPartonParticle',

    TagDir +'/LjetPtVsMinDRLjetTopPartonZoomDetector',
    TagDir +'/LjetPtVsMinDRLjetTopPartonZoomParticle',
    TagDir +'/LjetPtVsMinDRLjetTopPartonDetector',
    TagDir +'/LjetPtVsMinDRLjetTopPartonParticle',

    ]

names1dCuts = [

    [TagDir + '/MinDRLjetTopPartonParticle', TagDir + '/MinDRLjetTopPartonDetector'],
    [TagDir + '/MinDRLjetWPartonParticle', TagDir + '/MinDRLjetWPartonDetector'],
    [TagDir + '/MinDRLjetTopPartonZoomParticle', TagDir + '/MinDRLjetTopPartonZoomDetector'],
    [TagDir + '/MinDRLjetWPartonZoomParticle', TagDir + '/MinDRLjetWPartonZoomDetector'],
    
    [ LastCut + '/ParticleJetN', LastCut + '/DetectorJetN'],
    [ LastCut + '/ParticleJetPt', LastCut + '/DetectorJetPt'],
    [ LastCut + '/ParticleJet1Pt', LastCut + '/DetectorJet1Pt'],
    [ LastCut + '/ParticleJet2Pt', LastCut + '/DetectorJet2Pt'],
    [ LastCut + '/ParticleJet3Pt', LastCut + '/DetectorJet3Pt'],
    [ LastCut + '/ParticleJetRapidity', LastCut + '/DetectorJetRapidity'],
    [ LastCut + '/ParticleJet1Rapidity', LastCut + '/DetectorJet1Rapidity'],
    [ LastCut + '/ParticleJet2Rapidity', LastCut + '/DetectorJet2Rapidity'],
    [ LastCut + '/ParticleJet3Rapidity', LastCut + '/DetectorJet3Rapidity'],
    [ LastCut + '/ParticleJet4Rapidity', LastCut + '/DetectorJet4Rapidity'],
    [ LastCut + '/ParticleJetMass', LastCut + '/DetectorJetMass'],
    [ LastCut + '/ParticleJet1Mass', LastCut + '/DetectorJet1Mass'],
    [ LastCut + '/ParticleJet2Mass', LastCut + '/DetectorJet2Mass'],
    [ LastCut + '/ParticleJet3Mass', LastCut + '/DetectorJet3Mass'],
    [ LastCut + '/ParticleJet4Mass', LastCut + '/DetectorJet4Mass'],


    [ LastCut + '/ParticleLJetN', LastCut + '/DetectorLJetN'],
    [ LastCut + '/ParticleLJetPt', LastCut + '/DetectorLJetPt'],
    [ LastCut + '/ParticleLJet1Pt', LastCut + '/DetectorLJet1Pt'],
    [ LastCut + '/ParticleLJet2Pt', LastCut + '/DetectorLJet2Pt'],
    [ LastCut + '/ParticleLJet3Pt', LastCut + '/DetectorLJet3Pt'],
    [ LastCut + '/ParticleLJetRapidity', LastCut + '/DetectorLJetRapidity'],
    [ LastCut + '/ParticleLJet1Rapidity', LastCut + '/DetectorLJet1Rapidity'],
    [ LastCut + '/ParticleLJet2Rapidity', LastCut + '/DetectorLJet2Rapidity'],
    [ LastCut + '/ParticleLJet3Rapidity', LastCut + '/DetectorLJet3Rapidity'],
    [ LastCut + '/ParticleLJet4Rapidity', LastCut + '/DetectorLJet4Rapidity'],
    [ LastCut + '/ParticleLJetMass', LastCut + '/DetectorLJetMass'],
    [ LastCut + '/ParticleLJet1Mass', LastCut + '/DetectorLJet1Mass'],
    [ LastCut + '/ParticleLJet2Mass', LastCut + '/DetectorLJet2Mass'],
    [ LastCut + '/ParticleLJet3Mass', LastCut + '/DetectorLJet3Mass'],
    [ LastCut + '/ParticleLJet4Mass', LastCut + '/DetectorLJet4Mass'],
    [ LastCut + '/ParticleLJetTau32', LastCut + '/DetectorLJetTau32'],
    [ LastCut + '/ParticleLJet1Tau32', LastCut + '/DetectorLJet1Tau32'],
    [ LastCut + '/ParticleLJet2Tau32', LastCut + '/DetectorLJet2Tau32'],
    [ LastCut + '/ParticleLJet3Tau32', LastCut + '/DetectorLJet3Tau32'],
    [ LastCut + '/ParticleLJet4Tau32', LastCut + '/DetectorLJet4Tau32'],
    [ LastCut + '/ParticleLJetTau21', LastCut + '/DetectorLJetTau21'],
    [ LastCut + '/ParticleLJet1Tau21', LastCut + '/DetectorLJet1Tau21'],
    [ LastCut + '/ParticleLJet2Tau21', LastCut + '/DetectorLJet2Tau21'],
    [ LastCut + '/ParticleLJet3Tau21', LastCut + '/DetectorLJet3Tau21'],
    [ LastCut + '/ParticleLJet4Tau21', LastCut + '/DetectorLJet4Tau21'],

    
#    [ LastCut + '/ParticlebJetN', LastCut + '/DetectorbJetN'],
#    [ LastCut + '/ParticlebJetPt', LastCut + '/DetectorbJetPt'],
#    [ LastCut + '/ParticlebJetPt', LastCut + '/DetectorbJet1Pt'],
#    [ LastCut + '/ParticlebJetPt', LastCut + '/DetectorbJet2Pt'],
#    [ LastCut + '/ParticlebJetRapidity', LastCut + '/DetectorbJetRapidity'],
#    [ LastCut + '/ParticlebJetRapidity', LastCut + '/DetectorbJet1Rapidity'],
#    [ LastCut + '/ParticlebJetRapidity', LastCut + '/DetectorbJet2Rapidity'],


#    [ LastCut + '/ParticleLeptonPt', LastCut + '/DetectorLeptonPt'],
#    [ LastCut + '/ParticleNeutrinoPt', LastCut + '/DetectorNeutrinoPt'],


    ]


# controls:

names1dNoCuts = [
    [ FirstCut + '/ParticleJetN', FirstCut + '/DetectorJetN'],
    [ FirstCut + '/ParticleJetPt', FirstCut + '/DetectorJetPt'],
    [ FirstCut + '/ParticleJet1Pt', FirstCut + '/DetectorJet1Pt'],
    [ FirstCut + '/ParticleJet2Pt', FirstCut + '/DetectorJet2Pt'],
    [ FirstCut + '/ParticleJet3Pt', FirstCut + '/DetectorJet3Pt'],
    [ FirstCut + '/ParticleJetRapidity', FirstCut + '/DetectorJetRapidity'],
    [ FirstCut + '/ParticleJet1Rapidity', FirstCut + '/DetectorJet1Rapidity'],
    [ FirstCut + '/ParticleJet2Rapidity', FirstCut + '/DetectorJet2Rapidity'],
    [ FirstCut + '/ParticleJet3Rapidity', FirstCut + '/DetectorJet3Rapidity'],
    [ FirstCut + '/ParticleJet4Rapidity', FirstCut + '/DetectorJet4Rapidity'],
    [ FirstCut + '/ParticleJetMass', FirstCut + '/DetectorJetMass'],
    [ FirstCut + '/ParticleJet1Mass', FirstCut + '/DetectorJet1Mass'],
    [ FirstCut + '/ParticleJet2Mass', FirstCut + '/DetectorJet2Mass'],
    [ FirstCut + '/ParticleJet3Mass', FirstCut + '/DetectorJet3Mass'],
    [ FirstCut + '/ParticleJet4Mass', FirstCut + '/DetectorJet4Mass'],


    [ FirstCut + '/ParticleLJetN', FirstCut + '/DetectorLJetN'],
    [ FirstCut + '/ParticleLJetPt', FirstCut + '/DetectorLJetPt'],
    [ FirstCut + '/ParticleLJet1Pt', FirstCut + '/DetectorLJet1Pt'],
    [ FirstCut + '/ParticleLJet2Pt', FirstCut + '/DetectorLJet2Pt'],
    [ FirstCut + '/ParticleLJet3Pt', FirstCut + '/DetectorLJet3Pt'],
    [ FirstCut + '/ParticleLJetRapidity', FirstCut + '/DetectorLJetRapidity'],
    [ FirstCut + '/ParticleLJet1Rapidity', FirstCut + '/DetectorLJet1Rapidity'],
    [ FirstCut + '/ParticleLJet2Rapidity', FirstCut + '/DetectorLJet2Rapidity'],
    [ FirstCut + '/ParticleLJet3Rapidity', FirstCut + '/DetectorLJet3Rapidity'],
    [ FirstCut + '/ParticleLJet4Rapidity', FirstCut + '/DetectorLJet4Rapidity'],
    [ FirstCut + '/ParticleLJetMass', FirstCut + '/DetectorLJetMass'],
    [ FirstCut + '/ParticleLJet1Mass', FirstCut + '/DetectorLJet1Mass'],
    [ FirstCut + '/ParticleLJet2Mass', FirstCut + '/DetectorLJet2Mass'],
    [ FirstCut + '/ParticleLJet3Mass', FirstCut + '/DetectorLJet3Mass'],
    [ FirstCut + '/ParticleLJet4Mass', FirstCut + '/DetectorLJet4Mass'],
    [ FirstCut + '/ParticleLJetTau32', FirstCut + '/DetectorLJetTau32'],
    [ FirstCut + '/ParticleLJet1Tau32', FirstCut + '/DetectorLJet1Tau32'],
    [ FirstCut + '/ParticleLJet2Tau32', FirstCut + '/DetectorLJet2Tau32'],
    [ FirstCut + '/ParticleLJet3Tau32', FirstCut + '/DetectorLJet3Tau32'],
    [ FirstCut + '/ParticleLJet4Tau32', FirstCut + '/DetectorLJet4Tau32'],
    [ FirstCut + '/ParticleLJetTau21', FirstCut + '/DetectorLJetTau21'],
    [ FirstCut + '/ParticleLJet1Tau21', FirstCut + '/DetectorLJet1Tau21'],
    [ FirstCut + '/ParticleLJet2Tau21', FirstCut + '/DetectorLJet2Tau21'],
    [ FirstCut + '/ParticleLJet3Tau21', FirstCut + '/DetectorLJet3Tau21'],
    [ FirstCut + '/ParticleLJet4Tau21', FirstCut + '/DetectorLJet4Tau21'],

    
#    [ FirstCut + '/ParticlebJetN', FirstCut + '/DetectorbJetN'],
#    [ FirstCut + '/ParticlebJetPt', FirstCut + '/DetectorbJetPt'],
#    [ FirstCut + '/ParticlebJetPt', FirstCut + '/DetectorbJet1Pt'],
#    [ FirstCut + '/ParticlebJetPt', FirstCut + '/DetectorbJet2Pt'],
#    [ FirstCut + '/ParticlebJetRapidity', FirstCut + '/DetectorbJetRapidity'],
#    [ FirstCut + '/ParticlebJetRapidity', FirstCut + '/DetectorbJet1Rapidity'],
#    [ FirstCut + '/ParticlebJetRapidity', FirstCut + '/DetectorbJet2Rapidity'],


#    [ FirstCut + '/ParticleLeptonPt', FirstCut + '/DetectorLeptonPt'],
#    [ FirstCut + '/ParticleNeutrinoPt', FirstCut + '/DetectorNeutrinoPt'],


    ]

names2d_migras = [
# topology migrations:
    'SelectionMigrations',
# correlations:
    # old names before automation  
    # migrations:
    
    '2B0S/' + migradir + 'TopPt_Particle_Detector',
    '1B1S/' + migradir + 'TopPt_Particle_Detector',
    '0B2S/' + migradir + 'TopPt_Particle_Detector',
    'AnySel/' + migradir + 'TopPt_Particle_Detector',

    '2B0S/' + migradir + 'DiTopPt_Particle_Detector',
    '1B1S/' + migradir + 'DiTopPt_Particle_Detector',
    '0B2S/' + migradir + 'DiTopPt_Particle_Detector',
    'AnySel/' + migradir + 'DiTopPt_Particle_Detector',

    '2B0S/' + migradir + 'DiTopMass_Particle_Detector',
    '1B1S/' + migradir + 'DiTopMass_Particle_Detector',
    '0B2S/' + migradir + 'DiTopMass_Particle_Detector',
    'AnySel/' + migradir + 'DiTopMass_Particle_Detector',


    '2B0S/' + migradir + 'DiTopPout_Particle_Detector',
    '1B1S/' + migradir + 'DiTopPout_Particle_Detector',
    '0B2S/' + migradir + 'DiTopPout_Particle_Detector',
    'AnySel/' + migradir + 'DiTopPout_Particle_Detector',

    '2B0S/' + migradir + 'DiTopCosThetaStar_Particle_Detector',
    '1B1S/' + migradir + 'DiTopCosThetaStar_Particle_Detector',
    '0B2S/' + migradir + 'DiTopCosThetaStar_Particle_Detector',
    'AnySel/' + migradir + 'DiTopCosThetaStar_Particle_Detector',


    '2B0S/' + migradir + 'DiTopDeltaPhi_Particle_Detector',
    '1B1S/' + migradir + 'DiTopDeltaPhi_Particle_Detector',
    '0B2S/' + migradir + 'DiTopDeltaPhi_Particle_Detector',
    'AnySel/' + migradir + 'DiTopDeltaPhi_Particle_Detector',

    '2B0S/' + migradir + 'DiTopYboost_Particle_Detector',
    '1B1S/' + migradir + 'DiTopYboost_Particle_Detector',
    '0B2S/' + migradir + 'DiTopYboost_Particle_Detector',
    'AnySel/' + migradir + 'DiTopYboost_Particle_Detector',

    '2B0S/' + migradir + 'DiTopChittbar_Particle_Detector',
    '1B1S/' + migradir + 'DiTopChittbar_Particle_Detector',
    '0B2S/' + migradir + 'DiTopChittbar_Particle_Detector',
    'AnySel/' + migradir + 'DiTopChittbar_Particle_Detector',

    

    ]

#names1d = names1d_gen
# DEFAULT:
names1d = names1d_SBreco + names1d_Rel + names1d_Glob + names1dCuts + names1dNoCuts

# HACK!!!

#names1d = [    ['2B0S/ParticleTop1Pt' + suff,   '2B0S/DetectorTop1Pt' + suff],
#    ['2B0S/ParticleTop2Pt' + suff,   '2B0S/DetectorTop2Pt' + suff], 
#    ['2B0S/ParticleTopPt' + suff,    '2B0S/DetectorTopPt' + suff],
#]

#names1d = names1d_Glob[:1]
# HACK!!!
#names1d = [
#    ['1B1S/ParticleW1Mass',   '1B1S/DetectorW1Mass'],
#    ['1B1S/ParticleW2Mass',   '1B1S/DetectorW2Mass'], 
#    ['1B1S/ParticleWMass',    '1B1S/DetectorWMass'],
#    
#    ['0B2S/ParticleW1Mass',   '0B2S/DetectorW1Mass'],
#    ['0B2S/ParticleW2Mass',   '0B2S/DetectorW2Mass'], 
#    ['0B2S/ParticleWMass',    '0B2S/DetectorWMass'],
#]

# HACK!!!
#names2d = []
# names2d = names2d_all
# DEFAULT:
names2d = names2d_migras + names2d_tagging
#names2d = names2d_migras


