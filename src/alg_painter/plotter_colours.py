"""
This file is dedicated to the MERIAN unit colour palettes
developed by Karen van Neikerk and Charlotte Wright

I've converted most to hex_code rather than RGB tuples
"""

import matplotlib.pyplot as plt

MERIANBOW_4_MERIAN = {
    "M1": "#666666",
    "M2": "#710093",
    "M3": "#BC007B",
    "M4": "#EA005B",
    "M5": "#fc1d1d",
    "M6": "#FD9514",
    "M7": "#E9CB19",
    "M8": "#87BF13",
    "M9": "#00AB3E",
    "M10": "#00f2a1",
    "M11": "#005c66",
    "M12": "#00589E",
    "M13": "#006DDB",
    "M14": "#0080ff",
    "M15": "#A676FF",
    "M16": "#ff1aef",
    "M17": "#FF82CD",
    "M18": "#FF6B70",
    "M19": "#EE6A15",
    "M20": "#A16C00",
    "M21": "#4E6400",
    "M22": "#005200",
    "M23": "#00e251",
    "M24": "#009286",
    "M25": "#00B0E0",
    "M26": "#00C8FF",
    "M27": "#A1D0FF",
    "M28": "#ABA9E5",
    "M29": "#AE83B6",
    "M30": "#A56183",
    "M31": "#8E454F",
    "M32": "#6B3122",
    "MZ": "#404040FF",
}

MERIANBOW_MERIAN = {
    "M1": "#666666",
    "M2": "#D50062",
    "M3": "#F10059",
    "M4": "#FF104E",
    "M5": "#FF4F44",
    "M6": "#FF7839",
    "M7": "#FF9E2F",
    "M8": "#FF9700",
    "M9": "#CF8F00",
    "M10": "#9C8600",
    "M11": "#6A7A00",
    "M12": "#336C00",
    "M13": "#0A8500",
    "M14": "#009F00",
    "M15": "#00B846",
    "M16": "#00D27C",
    "M17": "#00EBB3",
    "M18": "#00D5BF",
    "M19": "#00BEC9",
    "M20": "#00A6D1",
    "M21": "#008FD5",
    "M22": "#0079D4",
    "M23": "#008BFC",
    "M24": "#009BFF",
    "M25": "#00A8FF",
    "M26": "#00B3FF",
    "M27": "#87BCFF",
    "M28": "#B198FF",
    "M29": "#C971FF",
    "M30": "#D546E1",
    "M31": "#D500B0",
    "M32": "#CC007E",
    "MZ": "#404040FF",
}

COLOUR_BLIND_MERIAN = {
    "M1": "#1F77B4",
    "M2": "#AEC7E8",
    "M3": "#FF7F0E",
    "M4": "#FFBB78",
    "M5": "#2CA02C",
    "M6": "#98DF8A",
    "M7": "#D62728",
    "M8": "#FF9896",
    "M9": "#9467BD",
    "M10": "#C5B0D5",
    "M11": "#8C564B",
    "M12": "#C49C94",
    "M13": "#E377C2",
    "M14": "#F7B6D2",
    "M15": "#7F7F7F",
    "M16": "#C7C7C7",
    "M17": "#008080",  # dark teal
    "M18": "#80BFBF",  # light teal
    "M19": "#17BECF",
    "M20": "#9EDAE5",
    "M21": "#393B79",
    "M22": "#5254A3",
    "M23": "#6B6ECF",
    "M24": "#9C9EDE",
    "M25": "#637939",
    "M26": "#8CA252",
    "M27": "#B5CF6B",
    "M28": "#CEDB9C",
    "M29": "#8C6D31",
    "M30": "#BD9E39",
    "M31": "#E7BA52",
    "MZ": "#404040FF",
}


def get_merian_colors_spectrum():
    """Generate colorblind-friendly spectrum using perceptually uniform colormaps"""
    turbo = plt.colormaps["turbo"]
    merian_elements = [f"M{i}" for i in range(1, 32)]
    colors = {}

    for i, merian in enumerate(merian_elements):
        position = 0.05 + (i / 30) * 0.90
        colors[merian] = turbo(position)

    colors["MZ"] = (0.25, 0.25, 0.25)
    return colors


MERIAN_COLOURS = {
    "spectrum": get_merian_colors_spectrum(),
    "merianbow": MERIANBOW_MERIAN,
    "merianbow4": MERIANBOW_4_MERIAN,
    "categorical": COLOUR_BLIND_MERIAN,
}
