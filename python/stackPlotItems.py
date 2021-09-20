#!/usr/bin/python

# jk 2020:
# for historical reason, histograms should be paired
# this is a remnant of shape cmp plotting when left was particle and right detector
# now each pair gets its own canvas
# so it does not matter what the pair is;)

LastCut = 'AnySel'
FirstCut = 'NoCuts'
TagDir = 'TagHistos'


ChiKeys = [ '2B0S',
            '1B1S',
            '0B2S',
            'AnySel']
FullKeys = [ 'NoCuts',
             '2B0S',
            '1B1S',
            '0B2S',
             'AnySel',]



# jk: note 29.5.2020:
# unf. no duplicities in histo names to be plotted are allowed as otherwiose TCanvas objects are deleted and the plotting code then crashes...

#########################################
def MakeGlobalPlotItems():
    items = []
    for topo in ChiKeys:
        items = items + [
            [topo + '/ParticleMet', topo + '/DetectorMet'],
            [topo + '/ParticleHTj', topo + '/DetectorHTj'],
            [topo + '/ParticleHTjPlusMet', topo + '/DetectorHTjPlusMet'],
            [topo + '/ParticleSumMj', topo + '/DetectorSumMj'],
            [topo + '/ParticleSumMJ', topo + '/DetectorSumMJ'],
            [topo + '/ParticleHTJ', topo + '/DetectorHTJ'],
            [topo + '/ParticlejAplanarity', topo + '/DetectorjAplanarity'],
            [topo + '/ParticlejSphericity', topo + '/DetectorjSphericity'],
            [topo + '/ParticleJAplanarity', topo + '/DetectorJAplanarity'],
            [topo + '/ParticleJSphericity', topo + '/DetectorJSphericity'],
            ]
    return items

#########################################
def MakeRelPlotItems():
    items = []
    for topo in ChiKeys:
        if topo == 'AnySel':
            continue
        items = items + [
            [topo + '/ParticleDiTopPtRel',   topo + '/DetectorDiTopPtRel'],
            [topo + '/ParticleDiTopPtGeo',   topo + '/DetectorDiTopPtGeo'],
            [topo + '/ParticleDiTopMassGeo',   topo + '/DetectorDiTopMassGeo'],
            [topo + '/ParticleTopPtRel',   topo + '/DetectorTopPtRel'],
            [topo + '/ParticleTop1PtRel',   topo + '/DetectorTop1PtRel'],
            [topo + '/ParticleTop2PtRel',   topo + '/DetectorTop2PtRel'],
            [topo + '/ParticleDiTopPoutRel',   topo + '/DetectorDiTopPoutRel'],
            [topo + '/ParticleDiTopPoutGeo',   topo + '/DetectorDiTopPoutGeo'],
        ]
    return items


#########################################
def MakeSpecialPlotItems(suff = '_denser'):
    items = []
    for topo in ChiKeys:
        if topo == 'AnySel':
            continue
        items = items + [
            [topo + '/ParticleDiTopPout' + suff,     topo + '/DetectorDiTopPout' + suff],
            [topo + '/ParticleDiTopCosThetaStar',    topo + '/DetectorDiTopCosThetaStar'],
            [topo + '/ParticleDiTopDelta',           topo + '/DetectorDiTopDelta'],
            [topo + '/ParticleDiTopYboost',          topo + '/DetectorDiTopYboost'],
            [topo + '/ParticleDiTopDeltaPhi',        topo + '/DetectorDiTopDeltaPhi'],
            [topo + '/ParticleDiTopChittbar',        topo + '/DetectorDiTopChittbar'],
            [topo + '/ParticleDiTopRttbar',          topo + '/DetectorDiTopRttbar'],
            
        ]
    return items

#########################################

def MakeNPlotItems():
    items = []
    for topo in ChiKeys:
        if topo == 'AnySel':
            continue
        items = items + [
            [topo + '/ParticleLJetN',   topo + '/DetectorLJetN'],
            [topo + '/ParticleJetN',   topo + '/DetectorJetN'],
            [topo + '/ParticlebJetN',   topo + '/DetectorbJetN'],
        ]
        if topo != 'NoCuts':
            items = items + [ [topo + '/ParticleTopN',   topo + '/DetectorTopN'],
                              [topo + '/ParticleWN',   topo + '/DetectorWN'],
            ]
    return items

#########################################

def MakeMassPlotItems(suff = '_denser'):
    items = []
    for topo in ChiKeys:
        if topo == 'AnySel':
            continue
        
        items = items + [
            [topo + '/ParticleLJetMass',   topo + '/DetectorLJetMass'],
            [topo + '/ParticleTop1Mass',   topo + '/DetectorTop1Mass'],
            [topo + '/ParticleTop2Mass',   topo + '/DetectorTop2Mass'], 
            [topo + '/ParticleTopMass',    topo + '/DetectorTopMass'],
            [topo + '/ParticleDiTopMass',    topo + '/DetectorDiTopMass'], #  + suff
            [topo + '/ParticleDiTopMass' + suff,    topo + '/DetectorDiTopMass' + suff], #  + suff
        ]


        if topo != '2B0S':
            items = items + [
                [topo + '/ParticleW1Mass',   topo + '/DetectorW1Mass'],
                [topo + '/ParticleWMass',    topo + '/DetectorWMass'],
                [topo + '/ParticleW1Pt',   topo + '/DetectorW1Pt'],
                [topo + '/ParticleWPt',    topo + '/DetectorWPt'],
                [topo + '/ParticleW1Rapidity',   topo + '/DetectorW1Rapidity'],
                [topo + '/ParticleWRapidity',    topo + '/DetectorWRapidity'],
                ]
        if topo == '0B2S' or topo == 'AnySel':
            items = items + [
                [ topo + '/ParticleW2Mass',   topo + '/DetectorW2Mass'],
                [ topo + '/ParticleW2Pt',   topo + '/DetectorW2Pt'],
                [ topo + '/ParticleW2Rapidity',   topo + '/DetectorW2Rapidity'], 
            ]
    return items

#########################################
def MakeMainKinemPlotItems(suff = '_denser'):
    items = []
    for topo in ChiKeys:
        if topo == 'AnySel':
            continue
        items = items + [
            [topo + '/ParticleDiTopPt' + suff,    topo + '/DetectorDiTopPt' + suff],
            [topo + '/ParticleTop1Pt' + suff,   topo + '/DetectorTop1Pt' + suff],
            [topo + '/ParticleTop2Pt' + suff,   topo + '/DetectorTop2Pt' + suff], 
            [topo + '/ParticleTopPt' + suff,    topo + '/DetectorTopPt' + suff],
            
            [topo + '/ParticleDiTopRapidity' + suff,    topo + '/DetectorDiTopRapidity' + suff],
            [topo + '/ParticleTop1Rapidity' + suff,   topo + '/DetectorTop1Rapidity' + suff],
            [topo + '/ParticleTop2Rapidity' + suff,   topo + '/DetectorTop2Rapidity' + suff], 
            [topo + '/ParticleTopRapidity' + suff,    topo + '/DetectorTopRapidity' + suff],
        ]
    return items

#########################################
def MakeCutPlotItems():
    items = []
    items = items + [
        [TagDir + '/MinDRLjetTopPartonParticle', TagDir + '/MinDRLjetTopPartonDetector'],
        [TagDir + '/MinDRLjetWPartonParticle', TagDir + '/MinDRLjetWPartonDetector'],
        [TagDir + '/MinDRLjetTopParton_zoomParticle', TagDir + '/MinDRLjetTopParton_zoomDetector'],
        [TagDir + '/MinDRLjetWParton_zoomParticle', TagDir + '/MinDRLjetWParton_zoomDetector'],
    ]

    items = items + [
        [ LastCut + '/ParticleJetN', LastCut + '/DetectorJetN'],
        [ LastCut + '/ParticleLJetN', LastCut + '/DetectorLJetN'],
        #    [ LastCut + '/ParticlebJetN', LastCut + '/DetectorbJetN'],
    ]

    ijets = [ '', '1', '2', '3', '4']
    for ijet in ijets:
        items = items + [
            [ LastCut + '/ParticleJet' + ijet + 'Pt', LastCut + '/DetectorJet' + ijet + 'Pt'],
            [ LastCut + '/ParticleJet' + ijet + 'Rapidity', LastCut + '/DetectorJet' + ijet + 'Rapidity'],
            [ LastCut + '/ParticleJet' + ijet + 'Mass', LastCut + '/DetectorJet' + ijet + 'Mass'],
            [ LastCut + '/ParticleLJet' + ijet + 'Pt', LastCut + '/DetectorLJet' + ijet + 'Pt'],
            [ LastCut + '/ParticleLJet' + ijet + 'Rapidity', LastCut + '/DetectorLJet' + ijet + 'Rapidity'],
            [ LastCut + '/ParticleLJet' + ijet + 'Mass', LastCut + '/DetectorLJet' + ijet + 'Mass'],
            [ LastCut + '/ParticleLJet' + ijet + 'Tau32', LastCut + '/DetectorLJet' + ijet + 'Tau32'],
            [ LastCut + '/ParticleLJet' + ijet + 'Tau21', LastCut + '/DetectorLJet' + ijet + 'Tau21'],
            [ LastCut + '/ParticlebJet' + ijet + 'Pt', LastCut + '/DetectorbJet' + ijet + 'Pt'],
            [ LastCut + '/ParticlebJet' + ijet + 'Rapidity', LastCut + '/DetectorbJet' + ijet + 'Rapidity'],
        ]

    return items

#########################################
def MakeFirstCutPlotItems():
    items = []
    # removed MinDR plots

    items = items + [
        [ FirstCut + '/ParticleLJetN', FirstCut + '/DetectorLJetN'],
        [ FirstCut + '/ParticlebJetN', FirstCut + '/DetectorbJetN'],
    ]

    ijets = [ '', '1', '2', '3', '4']
    for ijet in ijets:
        items = items + [
            [ FirstCut + '/ParticleJet' + ijet + 'Pt', FirstCut + '/DetectorJet' + ijet + 'Pt'],
            [ FirstCut + '/ParticleJet' + ijet + 'Rapidity', FirstCut + '/DetectorJet' + ijet + 'Rapidity'],
            [ FirstCut + '/ParticleJet' + ijet + 'Mass', FirstCut + '/DetectorJet' + ijet + 'Mass'],
            [ FirstCut + '/ParticleLJet' + ijet + 'Pt', FirstCut + '/DetectorLJet' + ijet + 'Pt'],
            [ FirstCut + '/ParticleLJet' + ijet + 'Rapidity', FirstCut + '/DetectorLJet' + ijet + 'Rapidity'],
            [ FirstCut + '/ParticleLJet' + ijet + 'Mass', FirstCut + '/DetectorLJet' + ijet + 'Mass'],
            [ FirstCut + '/ParticleLJet' + ijet + 'Tau32', FirstCut + '/DetectorLJet' + ijet + 'Tau32'],
            [ FirstCut + '/ParticleLJet' + ijet + 'Tau21', FirstCut + '/DetectorLJet' + ijet + 'Tau21'],
            [ FirstCut + '/ParticlebJet' + ijet + 'Pt', FirstCut + '/DetectorbJet' + ijet + 'Pt'],
            [ FirstCut + '/ParticlebJet' + ijet + 'Rapidity', FirstCut + '/DetectorbJet' + ijet + 'Rapidity'],
        ]
    return items

#########################################
