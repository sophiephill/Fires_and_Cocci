import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import geopandas as gpd


""" Returns number of fires in the dataframe """
def nunique(df):
    return len(df["hexid"].unique())


######## Cleaning functions #######################################################################


"""
Removes placebos / fires that had an MSE > cutoff. If a fire was fit poorly, all of its placebos are dropped.
Some fires are still outliers but have satisfactory pre MSE. Remove them using removeoutliers
"""

def removeBadfits(df, cutoff, removeoutliers = [], startMonth = 0):
    pre = df[df["month"] < startMonth]
    badfits = pre.groupby(["hexid","id"]).agg(lambda v: (v**2).mean())["gap"].reset_index()
    badfits = badfits[badfits["gap"] > cutoff].reset_index()
    for i in range(badfits.shape[0]):
        hexid = badfits.loc[i, "hexid"]
        ct = badfits.loc[i, "id"]
        ### If the badfit is a fire
        if "FR" in ct:
            df = df.drop(df[df["hexid"]==hexid].index)
        ### Else it's just a poorly fit placebo
        else:
            df = df.drop(df[(df["hexid"]==hexid) & (df["id"]==ct)].index)
    if removeoutliers:
        df = df.drop(df[df["hexid"].isin(removeoutliers)].index)
    return df

def smooth(df, col = "gap", by = 3, f = sum):
    smoothed = pd.Series([], dtype = "float64")
    for fire in df["hexid"].unique():
        sub = df[df["hexid"] == fire]
        if "id" in sub.columns:
            for ct in sub["id"].unique():
                sub2 = sub[sub["id"] == ct]
                gaps = sub2.rolling(by, on = "month", min_periods = 1).agg(f)[col]
                smoothed = pd.concat([smoothed, gaps], axis = 0)
        else:
            gaps = sub.rolling(by, on = "month", min_periods = 1).agg(f)[col]
            smoothed = pd.concat([smoothed, gaps], axis = 0)
    df2 = df.copy()
    df2[col] = smoothed
    return df2


def percentile(n):
    def percentile_(x):
        return np.nanpercentile(x, n)
    percentile_.__name__ = 'percentile_%s' % n
    return percentile_

"""
otherFires: stores fires that succeeded the fire so you can truncate plot effects plot before they overlapped
startMonth: for validation analysis when we varied the matching 
monthlyCi: View the bootstrapped ci for the treated effect
labelWHexid: if not, labels with fire name from output of NameFire
"""
def plotGsynGaps(df, size = "20km", otherfires = None, orangeFires = [], startMonth = 0, monthlyCi = None, putStarOn = []):
    
#     df["year"] = df["Fire.Label"].str.extract("(\d{4})")[0].astype(int)
    df = df.sort_values("igntn_d")
    
    nplots = len(df["hexid"].unique())
    nrows = -(-nplots // 4)
    fig, axs = plt.subplots(nrows, 4, sharex='col', sharey='row', figsize=(20,nrows*3),
                        gridspec_kw={'hspace': 0, 'wspace': 0})
    
    fig.suptitle("Month relative to ignition date", y = .1)
    index = 0
    rows = range(nrows)
    for hexid in df["hexid"].unique():
        sub = df.loc[df["hexid"]==hexid].sort_values(["hexid", "id", "month"])

        if nplots < 5:
            axis = axs[index]
        else:
            axis = axs[rows[index//4], index%4]
            
        pal = ["#b8b8b8" for s in sub["id"].unique()]
        
        ## Truncate post effects plot: end when another fire overlapped
        if (otherfires is not None) and (otherfires[otherfires["FID2"]==hexid].shape[0] > 0):
            whenoverlap = otherfires.loc[otherfires["FID2"]==hexid, "days_between"] // 30
            earliest = min(whenoverlap)
            sub = sub[sub["month"] <= earliest]    
        
        sns.lineplot(x = "month", y = "gap", data = sub, hue = "id", palette = pal, ax = axis,
                     legend = False, ci = None)
        col = "black"
        if hexid in orangeFires:    
            col = "#ff6600"
            
        sns.lineplot(x = "month", y = "gap", data = sub[sub["id"].str.contains("FR")], color = col,
                     ax = axis, legend = False, ci = None)
        
        if monthlyCi is not None:
            cis = monthlyCi.loc[monthlyCi["hexid"] == hexid, ["month", "CI.lower", "CI.upper"]]
            # If fires were truncated due to overlap in the post period, also truncate ribbon
            cis = cis[cis["month"].isin(sub["month"].unique())]
            axis.fill_between(cis["month"], cis["CI.lower"], cis["CI.upper"], color = c, alpha=0.2)      
        
        axis.axvline(startMonth, color = "black", linestyle = "dashed")
        axis.axhline(0, color = "black")
        axis.set_ylabel("")
        axis.set_xlabel("")
        if rows[index//4] == (nrows//2):
            axis.set_ylabel("Monthly differences in cases between the fire region and its synthetic control")

        label = sub["Fire.Label"].iloc[0]  
        if len(putStarOn) > 0 and (hexid in putStarOn):
            label = "*" + label + "*"
        if (index%4) == 0:
            yl = axis.get_ylim()[1] - 1
        axis.annotate(label, (min(sub["month"]), yl), weight = "bold")
        index += 1
        
"""
Takes in fire id of form FR45 and thexes dataframe and returns a name with the county, month, and year of fire formatted like Kern, Sept.2008
"""
def NameFire(ids, thexes):
    months = ["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sept", "oct", "nov", "dec"]
    names = []
    for fid in ids:
        if "ALL" in fid:
            names.append(fid)
            continue
        county = thexes.loc[thexes["HEXID"] == fid, "county"].str.replace("\sCounty", "").iloc[0]
        date = thexes.loc[thexes["HEXID"] == fid, "igntn_d"].iloc[0]
        m = int(date[5:7]) - 1
        name = county + ", " + months[m].title() + "." + date[:4]
        # Two fires the final 19 have the same name
        if fid == "FR41":
            name += " (2)"
        if fid == "FR7":
            name = 'San Joaquin, Jun.2009'
        names.append(name)
    return names


def getErrorbars(df):
    cimat = np.zeros(shape=(2, df.shape[0]))
    cimat[0] = abs(df["ATT"] - df["CI.lower"])
    cimat[1] = abs(df["ATT"] - df["CI.upper"])
    return cimat


def getPopDensity(population, returnValue = False):
    ## The 20km refers to the radius of the outer circle around the hexagon
    hexsize = 20 # 20,000 m
    sqrt3 = 3 ** 0.5
    R = hexsize * 2 / sqrt3
    hexareaSqKm = 1.5 * sqrt3 * (R ** 2)
    hexareaSqMiles = hexareaSqKm / 2.59
    density = population / hexareaSqMiles
    if returnValue:
        return round(density, 1)
    else:
        print( "{d:.0f} people per squared mile".format(d = density))




    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
