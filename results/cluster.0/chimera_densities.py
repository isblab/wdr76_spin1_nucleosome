import chimera
from chimera import openModels
from chimera import runCommand

# Names of proteins/domains for which we have created densities
prots = [
    "H2A0_n",
    "H2A0",
    "H2A0_c",
    "H2A1_n",
    "H2A1",
    "H2A1_c",
    "H2B0_n",
    "H2B0",
    "H2B1_n",
    "H2B1",
    "H30_n",
    "H30",
    "H31_n",
    "H31",
    "H40_n",
    "H40",
    "H41_n",
    "H41",
    "spin1_N",
    "spin1_tudor1",
    "spin1_td12",
    "spin1_tudor2",
    "spin1_td23",
    "spin1_tudor3",
    "wdr76_N",
    "wdr76_wd1",
    "wdr76_wd2",
    "wdr76_wd3",
    "wdr76_wd4",
    "wdr76_wd5",
    "wdr76_wd6",
    "wdr76_wd7",
]

# Set visualization thresholds
### About 20%
threshold = {
    "H2A0_n": 0.00666,
    "H2A0": 0.0634,
    "H2A0_c": 0.0106,
    "H2A1_n": 0.00912,
    "H2A1": 0.0656,
    "H2A1_c": 0.0081,
    "H2B0_n": 0.00866,
    "H2B0": 0.0682,
    "H2B1_n": 0.0106,
    "H2B1": 0.0756,
    "H30_n": 0.0108,
    "H30": 0.0874,
    "H31_n": 0.0284,
    "H31": 0.0646,
    "H40_n": 0.0256,
    "H40": 0.066,
    "H41_n": 0.00652,
    "H41": 0.0638,
    "spin1_N": 0.032,
    "spin1_tudor1": 0.068,
    "spin1_td12": 0.0268,
    "spin1_tudor2": 0.0384,
    "spin1_td23": 0.01342,
    "spin1_tudor3": 0.0316,
    "wdr76_N": 0.0302,
    "wdr76_wd1": 0.022,
    "wdr76_wd2": 0.0572,
    "wdr76_wd3": 0.0542,
    "wdr76_wd4": 0.0488,
    "wdr76_wd5": 0.03,
    "wdr76_wd6": 0.0336,
    "wdr76_wd7": 0.0306,
}


# Color of each protein/domain
# col_their = {
#     "H2A0":"#ffb300",
#     "H2A1":"#ffb300",
#     "H2B0":"#ffb300",
#     "H2B1":"#ffb300",
#     "H30":"#ffb300",
#     "H31":"#ffb300",
#     "H40":"#ffb300",
#     "H41":"#ffb300",
#     "spin1_N":"#e7fc47ae8a90",
#     "spin1_tudor1":"#e78a8a",
#     "spin1_tudor2":"#e7c18a",
#     "spin1_tudor3":"#e7e38a",
#     "spin1_C":"#e7fce3fb1402",
#     "wdr76_N":"#0000ff",
#     "wdr76_wd1":"#00000000ffff",
#     "wdr76_wd2":"#30c30000ffff",
#     "wdr76_wd3":"#616b0000ffff",
#     "wdr76_wd4":"#86180000ffff",
#     "wdr76_wd5":"#b6db0000ffff",
#     "wdr76_wd6":"#b6db0000cf3c",
#     "wdr76_wd7":"#b6db00005555",
#     "wdr76_C":"#ff0000",
# }


# ffb300
col = {
    "H2A0_n": "#fff75e",
    "H2A0": "gold",
    "H2A0_c": "#fdb833",
    "H2A1_n": "#fff75e",
    "H2A1": "gold",
    "H2A1_c": "#fdb833",
    "H2B0_n": "#f7a399",
    "H2B0": "#ef6351",
    "H2B1_n": "#f7a399",
    "H2B1": "#ef6351",
    "H30_n": "#ffc2d1",
    "H30": "#fb6f92",
    "H31_n": "#ffbfdf",
    "H31": "#ffd0bf",
    "H40_n": "#f4a57a",
    "H40": "#ba664b",
    "H41_n": "#f4a57a",
    "H41": "#ba664b",
    "spin1_N": "#96ED89",
    "spin1_tudor1": "#45BF55",
    "spin1_td12": "#45BF55",
    "spin1_tudor2": "#168039",
    "spin1_td23": "#168039",
    "spin1_tudor3": "#044D29",
    "wdr76_N": "#bcd9ea",
    "wdr76_wd1": "#8bbdd9",
    "wdr76_wd2": "#5ba4cf",
    "wdr76_wd3": "#298fca",
    "wdr76_wd4": "#0079bf",
    "wdr76_wd5": "#026aa7",
    "wdr76_wd6": "#055a8c",
    "wdr76_wd7": "#094c72",
}


runCommand("set bgcolor white")
i = 0

# Read localization density by component, both samples together
for p in prots:
    runCommand("open LPD_" + p + ".mrc")
    runCommand("volume #" + str(i) + " step 1 ")
    runCommand("volume #" + str(i) + " level " + str(threshold[p]))
    # runCommand('volume #'+str(i)+' transparency '+str(transp[p]))
    runCommand("color " + col[p] + " #" + str(i))
    # runCommand('2dlabels create ' + str(p) + '_lab text ' + str(p) + ' color '+col[p]+' size 30 xpos .1 ypos ' + str(0.9 - i / 25.0))
    i += 1

runCommand("open aligned_5gt0.pdb")
runCommand("open cluster_center_model.rmf3")
runCommand("color " + " dark magenta #" + str(i) + ":.I")
runCommand("color" + " dark magenta #" + str(i) + ":.J")

# runCommand("color " + " goldenrod #" + str(i) + ":.D")
# runCommand("color" + " goldenrod #" + str(i) + ":.H")

# runCommand("color " + " #fb6f92 #" + str(i) + ":.A")
# runCommand("color" + "  hot pink #" + str(i) + ":.E")

# runCommand("color " + " #ba664b #" + str(i) + ":.B")
# runCommand("color" + " #ba664b #" + str(i) + ":.F")

runCommand("~display #" + str(i))
runCommand("delete #" + str(i) + ":.A")
runCommand("delete #" + str(i) + ":.B")
runCommand("delete #" + str(i) + ":.C")
runCommand("delete #" + str(i) + ":.D")
runCommand("delete #" + str(i) + ":.E")
runCommand("delete #" + str(i) + ":.F")
runCommand("delete #" + str(i) + ":.G")
runCommand("delete #" + str(i) + ":.H")


"""
H2A: #ff935c #c8b6a6
H2B: #ff6d3f #8d7b68
H3.0: #ffadbc
H3.1: #ec7272
H4: #e5ba73 #faab78
SPIN1: #a9eca2
"""
