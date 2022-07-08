#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 15:03:07 2021

@author: Dartoon
"""

ID_list = [
    '022404.85+014941.9',  #!!! pick it
    '022906.04-051428.9',  #!!! pick it
    '092532.13-020806.1',  #!!! pick it, PS0 fainter but similar color
    '095218.04-000459.1', #SAVED
    '105458.01+043310.6',  #!!! pick it, PS0 slightly redder. PS1 too faint
    '122144.31-004144.1', #SAVED
    '124618.51-001750.2',  #!!! pick it, very similar and bright color
    '150216.66+025719.8',         #No R band data but looks interesting.
    '162501.98+430931.6', #!!! Pick it 
    '220642.82+003016.2', #!!! pick it, close and color similar.
    '230402.77-003855.4', #Very close, PS1 redder
    # '090559.70+022408.2', #SAVED However, too faint
    # '141911.59+520545.2', #!!! John's proposal, But missed I band in S20a, but have I band in S19a  
    # '013834.18-000509.3',  #similar color PS1 fainter
    # '090654.53+021315.2',  #!!! Gaia detect parallax, PS0 redder. Probably not select it  #!!! In John's propsal
    # # Close to star line:
    # '144308.16-004913.4', #!!! pick it, PS1 slight red.  Shenli is running on Gemini
    # '145347.46+003927.0', #!!! pick it , Shenli is running on Gemini
    # '001459.72+002319.2', #interesting, but PS1 too close to star line.
    # '013736.57+005742.3',  #G band almost no detection for PS1,  ##### Too close to star line
    # # Too red or faint:
    # '023829.90-011224.2',  #Dual feature not clear. Too close? Also very faint.
    # '104644.31+000329.7',     #Too faint for both. but Similr color
    # '134257.16-013912.9',        #G band no detection, but PS0 also red. 
    # '220422.46+073138.3', #Red shift very high, PS1 only detected in I band
    # '011227.87-003151.6',     #G band  almost no detection
    # '012110.93+010703.3',     #G band no detection
    # '023600.28-010432.3',       #G band almost no detection, PS1 too tiny.
    # '131512.46+015021.6',     #G band no detection, but GAIA get nan.
    # '221011.62-001654.9',     #PS1 red and G band no detection.
    # '132441.58-015401.8',      #PS1 too red
    # '152112.96+441452.5',   #PS1 too red  
    # '153008.91+425634.8',         #Faint, PS1 G band no detection.
    # '233718.07+002550.6',     #PS1 too red
    # '120417.10+003653.7',     #PS1 too red, PS1 G almost no detection 
    # '124604.03-010954.6',     #PS1 too red,
    # '125216.06+003141.1',     #PS1 too red,    #In John's proposal
    # '221101.45+001449.0',     #PS1 too red and PS1 host can be seen. PS1 also faint in G band.  #In John's proposal
    # '220910.38-001601.5',     #PS1 to red, almost no detection.   #In John's proposal
    # '011935.29-002033.5',  #John's proposal, z~0.86, But too red.
    # '144034.78+441520.5', #John's proposal, z~0.8, But too red.
    #Observed error：
    #####'095218.04-000459.1',  # pick it, similar color, Only in I band. Observed error?
    ##### '000050.56-013055.2',  #Probably a error on I band, see also a similar feature on the bottom star.
    ##### '003659.26-001922.7',   ##Only in I band, They are measure error, see also a similar feature on the bottom star.
    ##### '003659.43-001850.2',     ##Only in I band, see also a similar feature on the bottom star.
  ]
# ID_list.sort()

#All the ones with z>1.3 and sep < 0.7"
# ID_list = ['022404.85+014941.9',
#  '022906.04-051428.9',
#  '092532.13-020806.1',
#  '105458.01+043310.6',
#  '124618.51-001750.2',
#  '150216.66+025719.8',
#  '162501.98+430931.6',
#  '220642.82+003016.2',
#  '230402.77-003855.4',
#  '022509.04-001346.4',
#  '013930.80-002131.6',
#  '010553.76+000232.1',
#  '023855.60-000750.1',
#  '013526.15+004135.8',
#  '003828.78+001000.1',
#  '095218.04-000459.1',
#  '122144.31-004144.1',
#  '003659.43-001850.2',
#  '015741.34+023911.8',
#  '003659.26-001922.7',
#  '002938.47+002039.7',
#  '152445.62+440949.5',
#  '023829.90-011224.2',
#  '013327.34-003612.7',
#  '090559.70+022408.2']
