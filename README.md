
## Pymaceuticals, Inc.

#### OBSERVED TRENDS:

1. Between the three treatment drugs (Capomulin, Infubinol, Ketapril) and Placebo; Ketapril has the worst effect in tumor volume growth over the course of a 45-day clinical trial period.  
<br>
2. Capomulin, has a noticeable reduction effect on the tumor volume over a span of 45 days. The spread of metastatic sites is also less for mice given Capomulin and there are significantly fewer tumors detected.  
<br>
3. According to the __Survival During Treatment__, slightly more mice given the drug Infubinol died in comparison to those given Ketapril. On the other hand, the reports show that mice given Ketapril had larger tumor growth and metastatic site counts. These results suggest that tumor size and spread alone is not a determining factor in the count of mouse deaths.


```python
# dependencies
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
```

I used sem() in error bars [`pandas.DataFrame.sem`](http://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.sem.html)


```python
# requires scipy.stats import sem
from scipy.stats import sem
```

#### Drug Trial Analysis


```python
# read csv files
df_ct_data = pd.read_csv("raw_data/clinicaltrial_data.csv")
df_md_data = pd.read_csv("raw_data/mouse_drug_data.csv")

# merge datasets on Mouse ID to get drug per trial
df_ct_drug_data = pd.merge(df_ct_data, df_md_data, on="Mouse ID")
df_ct_drug_data = df_ct_drug_data.sort_values(["Timepoint", "Tumor Volume (mm3)"])

# reset index to remove original index column, change data view
df_ct_drug_data = df_ct_drug_data.reset_index()
df_ct_drug_data.drop(df_ct_drug_data.columns[0], axis=1, inplace=True)
df_ct_drug_data.head(10)
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Mouse ID</th>
      <th>Timepoint</th>
      <th>Tumor Volume (mm3)</th>
      <th>Metastatic Sites</th>
      <th>Drug</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>b128</td>
      <td>0</td>
      <td>45.0</td>
      <td>0</td>
      <td>Capomulin</td>
    </tr>
    <tr>
      <th>1</th>
      <td>f932</td>
      <td>0</td>
      <td>45.0</td>
      <td>0</td>
      <td>Ketapril</td>
    </tr>
    <tr>
      <th>2</th>
      <td>g107</td>
      <td>0</td>
      <td>45.0</td>
      <td>0</td>
      <td>Ketapril</td>
    </tr>
    <tr>
      <th>3</th>
      <td>a457</td>
      <td>0</td>
      <td>45.0</td>
      <td>0</td>
      <td>Ketapril</td>
    </tr>
    <tr>
      <th>4</th>
      <td>c819</td>
      <td>0</td>
      <td>45.0</td>
      <td>0</td>
      <td>Ketapril</td>
    </tr>
    <tr>
      <th>5</th>
      <td>h246</td>
      <td>0</td>
      <td>45.0</td>
      <td>0</td>
      <td>Ketapril</td>
    </tr>
    <tr>
      <th>6</th>
      <td>p189</td>
      <td>0</td>
      <td>45.0</td>
      <td>0</td>
      <td>Ketapril</td>
    </tr>
    <tr>
      <th>7</th>
      <td>n923</td>
      <td>0</td>
      <td>45.0</td>
      <td>0</td>
      <td>Ketapril</td>
    </tr>
    <tr>
      <th>8</th>
      <td>q119</td>
      <td>0</td>
      <td>45.0</td>
      <td>0</td>
      <td>Ketapril</td>
    </tr>
    <tr>
      <th>9</th>
      <td>f993</td>
      <td>0</td>
      <td>45.0</td>
      <td>0</td>
      <td>Naftisol</td>
    </tr>
  </tbody>
</table>
</div>




```python
# define variables for drug names/column headers and scatter plot colors
# Capomulin, Infubinol, Ketapril, and Placebo
dr1 = "Capomulin"
dr2 = "Infubinol"
dr3 = "Ketapril"
placebo = "Placebo"

# colors
color_rd = "red"
color_bl = "blue"
color_gr = "green"
color_blk = "black"
```


```python
# DataFrames to hold data samples for standard error calculations 
# measurements for each drug at each timepoint or nested lists
df_capomulin_data = df_ct_drug_data[df_ct_drug_data["Drug"] == "Capomulin"]
df_infubinol_data = df_ct_drug_data[df_ct_drug_data["Drug"] == "Infubinol"]
df_ketapril_data = df_ct_drug_data[df_ct_drug_data["Drug"] == "Ketapril"]
df_placebo_data = df_ct_drug_data[df_ct_drug_data["Drug"] == "Placebo"]

# list to hold timepoint values
list_timepts = df_ct_drug_data["Timepoint"].unique()

# Tumor Volume data sets for timepoints
tv_dr1 = [df_capomulin_data.loc[df_capomulin_data["Timepoint"] == tp, "Tumor Volume (mm3)"] for tp in list_timepts]
tv_dr2 = [df_infubinol_data.loc[df_infubinol_data["Timepoint"] == tp, "Tumor Volume (mm3)"] for tp in list_timepts]
tv_dr3 = [df_ketapril_data.loc[df_ketapril_data["Timepoint"] == tp, "Tumor Volume (mm3)"] for tp in list_timepts]
tv_placebo = [df_placebo_data.loc[df_placebo_data["Timepoint"] == tp, "Tumor Volume (mm3)"] for tp in list_timepts]

# Metastatic Sites data sets for timepoints 
ms_dr1 = [df_capomulin_data.loc[df_capomulin_data["Timepoint"] == tp, "Metastatic Sites"] for tp in list_timepts]
ms_dr2 = [df_infubinol_data.loc[df_infubinol_data["Timepoint"] == tp, "Metastatic Sites"] for tp in list_timepts]
ms_dr3 = [df_ketapril_data.loc[df_ketapril_data["Timepoint"] == tp, "Metastatic Sites"] for tp in list_timepts]
ms_placebo = [df_placebo_data.loc[df_placebo_data["Timepoint"] == tp, "Metastatic Sites"] for tp in list_timepts]
```

#### Custom Helper Functions:

__drug_summary_table__
* Outputs summary table for our 3 drugs and placebo
* DataFrame of average values by timepoint, column to measure

__scatter_plots__
* DataFrame of data summary, measurement, plot title, drug names, and placebo
* Determine color of markers based on measured tumor growth or reduction 
   


```python
def drug_summary_table(df, col):

    # sort by drugs we care about
    df_dr1 = df[df["Drug"] == "Capomulin"]
    df_dr2 = df[df["Drug"] == "Infubinol"]
    df_dr3 = df[df["Drug"] == "Ketapril"]
    df_placebo = df[df["Drug"] == "Placebo"]
        
    # merge and rename
    df_drug_summary = pd.merge(df_dr1, df_dr2, on="Timepoint")
    df_drug_summary = df_drug_summary.rename(columns={
                                                    f"{col}_x": "Capomulin",
                                                    f"{col}_y": "Infubinol"
                                                    }
                                            )
    df_drug_summary = pd.merge(df_drug_summary, df_dr3, on="Timepoint")
    df_drug_summary = pd.merge(df_drug_summary, df_placebo, on="Timepoint")
    df_drug_summary = df_drug_summary.rename(columns={
                                                    f"{col}_x": "Ketapril",
                                                    f"{col}_y": "Placebo"
                                                    }
                                            )
    df_drug_summary.drop(df_drug_summary.columns[[0,3,5,7]], axis=1, inplace=True)
    return df_drug_summary
```


```python
def scatter_plots(df, m, plot_title):
    
    # plot series data for timepoints, including each drug's Tumor or Metastatic Site changes over time
    timepoints = df.loc[:, "Timepoint"]
    s_dr1 = df.loc[:, dr1]
    s_dr2 = df.loc[:, dr2]
    s_dr3  = df.loc[:, dr3]
    s_placebo = df.loc[:, placebo]
    
    # set ticks, limits, labels, and title
    fig, ax = plt.subplots()
    ax.set_xticks(np.arange(min(timepoints), max(timepoints)+1, 5))
    ax.tick_params(direction="out", color="black", width=1, length=5, axis="both", pad=2)
    ax.set_xlim(-1, max(timepoints) + 3)
    ax.set_xlabel("Time (Days)", fontsize=12)
    ax.set_ylabel(m, fontsize=12)   
    ax.set_title(plot_title, fontsize=18)
    ax.set_xmargin = 20
    
    # add padding below and to left of the ticks in the axes
    ax.xaxis.labelpad = 10
    ax.yaxis.labelpad = 10
    ax.title.set_position([.5, 1.04])
    
    # scatter plots, assign handles for our legend based on colors above
    handle1 = ax.scatter(timepoints, s_dr1, facecolors=color_rd, marker="o", label=dr1)
    handle2 = ax.scatter(timepoints, s_dr2, facecolors=color_bl, marker="^", label=dr2)
    handle3 = ax.scatter(timepoints, s_dr3, facecolors=color_gr,  marker="s", label=dr3)
    handle4 = ax.scatter(timepoints, s_placebo, facecolors=color_blk, marker="D", label=placebo)
    
    # Build line plots on top of scatters with same colors as above
    ax.plot(timepoints, s_dr1, '--', color=color_rd, linewidth=1)
    ax.plot(timepoints, s_dr2, '--', color=color_bl, linewidth=1)
    ax.plot(timepoints, s_dr3, '--', color=color_gr, linewidth=1)
    ax.plot(timepoints, s_placebo, '--', color=color_blk, linewidth=1)

    # define scatter plot legend
    ax.legend(handles=[handle1, handle2, handle3, handle4], loc='best', frameon=True, facecolor="white")

    # Build error bar means and standard error datasets + Define means and standard error margin values
    if (m != "Mouse Count"):
        
        if (m == "Tumor Volume (mm3)"):
            means_dr1   = [np.mean(c) for c in tv_dr1]
            means_dr2   = [np.mean(i) for i in tv_dr2]
            means_dr3   = [np.mean(k) for k in tv_dr3]
            means_placebo = [np.mean(p) for p in tv_placebo]
            se_dr1      = [sem(c) for c in tv_dr1]
            se_dr2      = [sem(i) for i in tv_dr2]
            se_dr3      = [sem(k) for k in tv_dr3]
            se_placebo    = [sem(p) for p in tv_placebo]
        elif (m == "Metastatic Sites"):
            means_dr1   = [np.mean(c) for c in ms_dr1]
            means_dr2   = [np.mean(i) for i in ms_dr2]
            means_dr3   = [np.mean(k) for k in ms_dr3]
            means_placebo = [np.mean(p) for p in ms_placebo]
            se_dr1      = [sem(c) for c in ms_dr1]
            se_dr2      = [sem(i) for i in ms_dr2]
            se_dr3      = [sem(k) for k in ms_dr3]
            se_placebo    = [sem(p) for p in ms_placebo]

        # Place error bars in axis
        ax.errorbar(timepoints, means_dr1, se_dr1, color=color_rd, fmt="o", capsize=5, elinewidth=1, markeredgewidth=1)
        ax.errorbar(timepoints, means_dr2, se_dr2, color=color_bl, fmt="^", capsize=5, elinewidth=1, markeredgewidth=1)
        ax.errorbar(timepoints, means_dr3, se_dr3, color=color_gr, fmt="s", capsize=5, elinewidth=1, markeredgewidth=1)
        ax.errorbar(timepoints, means_placebo, se_placebo, color=color_blk, fmt="D", capsize=5, elinewidth=1, markeredgewidth=1)
```

#### Table: Tumor Volume Relation To Drug Treatment


```python
# set variable
column_to_measure = "Tumor Volume (mm3)"

# group dataframe to display average tumor growth over the course of the timepoints
drugs_grouping = df_ct_drug_data.groupby(["Drug", "Timepoint"])

# List the average Tumor Volume for groups in dataframe
df_avg_tumorvolume_by_drug = pd.DataFrame(drugs_grouping[column_to_measure].mean()).reset_index()
df_avg_tumorvolume_by_drug.head(10)
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Drug</th>
      <th>Timepoint</th>
      <th>Tumor Volume (mm3)</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Capomulin</td>
      <td>0</td>
      <td>45.000000</td>
    </tr>
    <tr>
      <th>1</th>
      <td>Capomulin</td>
      <td>5</td>
      <td>44.266086</td>
    </tr>
    <tr>
      <th>2</th>
      <td>Capomulin</td>
      <td>10</td>
      <td>43.084291</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Capomulin</td>
      <td>15</td>
      <td>42.064317</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Capomulin</td>
      <td>20</td>
      <td>40.716325</td>
    </tr>
    <tr>
      <th>5</th>
      <td>Capomulin</td>
      <td>25</td>
      <td>39.939528</td>
    </tr>
    <tr>
      <th>6</th>
      <td>Capomulin</td>
      <td>30</td>
      <td>38.769339</td>
    </tr>
    <tr>
      <th>7</th>
      <td>Capomulin</td>
      <td>35</td>
      <td>37.816839</td>
    </tr>
    <tr>
      <th>8</th>
      <td>Capomulin</td>
      <td>40</td>
      <td>36.958001</td>
    </tr>
    <tr>
      <th>9</th>
      <td>Capomulin</td>
      <td>45</td>
      <td>36.236114</td>
    </tr>
  </tbody>
</table>
</div>



#### Summary Table


```python
# Create drug summary table for Tumor Volume
df_volume_summary = drug_summary_table(df_avg_tumorvolume_by_drug, column_to_measure)
df_volume_summary.head(10)
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Timepoint</th>
      <th>Capomulin</th>
      <th>Infubinol</th>
      <th>Ketapril</th>
      <th>Placebo</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0</td>
      <td>45.000000</td>
      <td>45.000000</td>
      <td>45.000000</td>
      <td>45.000000</td>
    </tr>
    <tr>
      <th>1</th>
      <td>5</td>
      <td>44.266086</td>
      <td>47.062001</td>
      <td>47.389175</td>
      <td>47.125589</td>
    </tr>
    <tr>
      <th>2</th>
      <td>10</td>
      <td>43.084291</td>
      <td>49.403909</td>
      <td>49.582269</td>
      <td>49.423329</td>
    </tr>
    <tr>
      <th>3</th>
      <td>15</td>
      <td>42.064317</td>
      <td>51.296397</td>
      <td>52.399974</td>
      <td>51.359742</td>
    </tr>
    <tr>
      <th>4</th>
      <td>20</td>
      <td>40.716325</td>
      <td>53.197691</td>
      <td>54.920935</td>
      <td>54.364417</td>
    </tr>
    <tr>
      <th>5</th>
      <td>25</td>
      <td>39.939528</td>
      <td>55.715252</td>
      <td>57.678982</td>
      <td>57.482574</td>
    </tr>
    <tr>
      <th>6</th>
      <td>30</td>
      <td>38.769339</td>
      <td>58.299397</td>
      <td>60.994507</td>
      <td>59.809063</td>
    </tr>
    <tr>
      <th>7</th>
      <td>35</td>
      <td>37.816839</td>
      <td>60.742461</td>
      <td>63.371686</td>
      <td>62.420615</td>
    </tr>
    <tr>
      <th>8</th>
      <td>40</td>
      <td>36.958001</td>
      <td>63.162824</td>
      <td>66.068580</td>
      <td>65.052675</td>
    </tr>
    <tr>
      <th>9</th>
      <td>45</td>
      <td>36.236114</td>
      <td>65.755562</td>
      <td>70.662958</td>
      <td>68.084082</td>
    </tr>
  </tbody>
</table>
</div>



#### Scatter Plot


```python
# Call function to create scatter plots for Tumor Volume data and disply data
scatter_plots(df_volume_summary, column_to_measure, "Tumor Response to Treatment")
plt.show()
```


![png](Images/output_17_0.png)


#### Table: Average Metastatic Sites By Timepoint & Response To Treatment


```python
column_to_measure = "Metastatic Sites"

df_avg_mssites_by_drug = pd.DataFrame(drugs_grouping[column_to_measure].mean()).reset_index()
df_avg_mssites_by_drug.head(10)
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Drug</th>
      <th>Timepoint</th>
      <th>Metastatic Sites</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Capomulin</td>
      <td>0</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>1</th>
      <td>Capomulin</td>
      <td>5</td>
      <td>0.160000</td>
    </tr>
    <tr>
      <th>2</th>
      <td>Capomulin</td>
      <td>10</td>
      <td>0.320000</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Capomulin</td>
      <td>15</td>
      <td>0.375000</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Capomulin</td>
      <td>20</td>
      <td>0.652174</td>
    </tr>
    <tr>
      <th>5</th>
      <td>Capomulin</td>
      <td>25</td>
      <td>0.818182</td>
    </tr>
    <tr>
      <th>6</th>
      <td>Capomulin</td>
      <td>30</td>
      <td>1.090909</td>
    </tr>
    <tr>
      <th>7</th>
      <td>Capomulin</td>
      <td>35</td>
      <td>1.181818</td>
    </tr>
    <tr>
      <th>8</th>
      <td>Capomulin</td>
      <td>40</td>
      <td>1.380952</td>
    </tr>
    <tr>
      <th>9</th>
      <td>Capomulin</td>
      <td>45</td>
      <td>1.476190</td>
    </tr>
  </tbody>
</table>
</div>



#### Summary Table of Metastatic Sites by Timepoint (w/ columns for treatment drugs)


```python
df_mssites_summary = drug_summary_table(df_avg_mssites_by_drug, column_to_measure)
df_mssites_summary.head(10)
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Timepoint</th>
      <th>Capomulin</th>
      <th>Infubinol</th>
      <th>Ketapril</th>
      <th>Placebo</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>1</th>
      <td>5</td>
      <td>0.160000</td>
      <td>0.280000</td>
      <td>0.304348</td>
      <td>0.375000</td>
    </tr>
    <tr>
      <th>2</th>
      <td>10</td>
      <td>0.320000</td>
      <td>0.666667</td>
      <td>0.590909</td>
      <td>0.833333</td>
    </tr>
    <tr>
      <th>3</th>
      <td>15</td>
      <td>0.375000</td>
      <td>0.904762</td>
      <td>0.842105</td>
      <td>1.250000</td>
    </tr>
    <tr>
      <th>4</th>
      <td>20</td>
      <td>0.652174</td>
      <td>1.050000</td>
      <td>1.210526</td>
      <td>1.526316</td>
    </tr>
    <tr>
      <th>5</th>
      <td>25</td>
      <td>0.818182</td>
      <td>1.277778</td>
      <td>1.631579</td>
      <td>1.941176</td>
    </tr>
    <tr>
      <th>6</th>
      <td>30</td>
      <td>1.090909</td>
      <td>1.588235</td>
      <td>2.055556</td>
      <td>2.266667</td>
    </tr>
    <tr>
      <th>7</th>
      <td>35</td>
      <td>1.181818</td>
      <td>1.666667</td>
      <td>2.294118</td>
      <td>2.642857</td>
    </tr>
    <tr>
      <th>8</th>
      <td>40</td>
      <td>1.380952</td>
      <td>2.100000</td>
      <td>2.733333</td>
      <td>3.166667</td>
    </tr>
    <tr>
      <th>9</th>
      <td>45</td>
      <td>1.476190</td>
      <td>2.111111</td>
      <td>3.363636</td>
      <td>3.272727</td>
    </tr>
  </tbody>
</table>
</div>



#### Scatter Plot: Tumor Response to Treatment Over Time


```python
scatter_plots(df_mssites_summary, column_to_measure, "Metastatic Spread During Treatment")
plt.show()
```


![png](Images/output_23_0.png)


#### Survival Rate: Table Averaging Mouse Count by Timepoint


```python
# repeat tables for scatter plot
df_avg_mice_by_drug = pd.DataFrame(drugs_grouping["Tumor Volume (mm3)"].count()).reset_index()
df_avg_mice_by_drug = df_avg_mice_by_drug.rename(columns={f"Tumor Volume (mm3)": "Mouse Count"})

df_avg_mice_by_drug.head(10)
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Drug</th>
      <th>Timepoint</th>
      <th>Mouse Count</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Capomulin</td>
      <td>0</td>
      <td>25</td>
    </tr>
    <tr>
      <th>1</th>
      <td>Capomulin</td>
      <td>5</td>
      <td>25</td>
    </tr>
    <tr>
      <th>2</th>
      <td>Capomulin</td>
      <td>10</td>
      <td>25</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Capomulin</td>
      <td>15</td>
      <td>24</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Capomulin</td>
      <td>20</td>
      <td>23</td>
    </tr>
    <tr>
      <th>5</th>
      <td>Capomulin</td>
      <td>25</td>
      <td>22</td>
    </tr>
    <tr>
      <th>6</th>
      <td>Capomulin</td>
      <td>30</td>
      <td>22</td>
    </tr>
    <tr>
      <th>7</th>
      <td>Capomulin</td>
      <td>35</td>
      <td>22</td>
    </tr>
    <tr>
      <th>8</th>
      <td>Capomulin</td>
      <td>40</td>
      <td>21</td>
    </tr>
    <tr>
      <th>9</th>
      <td>Capomulin</td>
      <td>45</td>
      <td>21</td>
    </tr>
  </tbody>
</table>
</div>



#### Summary Table: Average Mouse Count by Timepoint (w/ columns for treatment drugs)


```python
df_mscount_summary = drug_summary_table(df_avg_mice_by_drug, "Mouse Count")
df_mscount_summary.head(10)
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Timepoint</th>
      <th>Capomulin</th>
      <th>Infubinol</th>
      <th>Ketapril</th>
      <th>Placebo</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0</td>
      <td>25</td>
      <td>25</td>
      <td>25</td>
      <td>25</td>
    </tr>
    <tr>
      <th>1</th>
      <td>5</td>
      <td>25</td>
      <td>25</td>
      <td>23</td>
      <td>24</td>
    </tr>
    <tr>
      <th>2</th>
      <td>10</td>
      <td>25</td>
      <td>21</td>
      <td>22</td>
      <td>24</td>
    </tr>
    <tr>
      <th>3</th>
      <td>15</td>
      <td>24</td>
      <td>21</td>
      <td>19</td>
      <td>20</td>
    </tr>
    <tr>
      <th>4</th>
      <td>20</td>
      <td>23</td>
      <td>20</td>
      <td>19</td>
      <td>19</td>
    </tr>
    <tr>
      <th>5</th>
      <td>25</td>
      <td>22</td>
      <td>18</td>
      <td>19</td>
      <td>17</td>
    </tr>
    <tr>
      <th>6</th>
      <td>30</td>
      <td>22</td>
      <td>17</td>
      <td>18</td>
      <td>15</td>
    </tr>
    <tr>
      <th>7</th>
      <td>35</td>
      <td>22</td>
      <td>12</td>
      <td>17</td>
      <td>14</td>
    </tr>
    <tr>
      <th>8</th>
      <td>40</td>
      <td>21</td>
      <td>10</td>
      <td>15</td>
      <td>12</td>
    </tr>
    <tr>
      <th>9</th>
      <td>45</td>
      <td>21</td>
      <td>9</td>
      <td>11</td>
      <td>11</td>
    </tr>
  </tbody>
</table>
</div>



#### Scatter Plot: Mouse Count Over Time for Treatment Drugs


```python
# use same scatter_plot function
scatter_plots(df_mscount_summary, "Mouse Count", "Survival During Treatment")
plt.show()
```


![png](Images/output_29_0.png)


#### Bar Graph: Comparing Total % Tumor Volume Change for Treatment Drugs (in 45 days)


```python
dr1_percent = (df_volume_summary[dr1][len(list_timepts)-1] - df_volume_summary[dr1][0]) / df_volume_summary[dr1][0]
dr2_percent = (df_volume_summary[dr2][len(list_timepts)-1] - df_volume_summary[dr2][0]) / df_volume_summary[dr2][0]
dr3_percent  = (df_volume_summary[dr3][len(list_timepts)-1] - df_volume_summary[dr3][0]) / df_volume_summary[dr3][0]
placebo_percent   = (df_volume_summary[placebo][len(list_timepts)-1] - df_volume_summary[placebo][0])  / df_volume_summary[placebo][0]

percent_change = pd.Series([dr1_percent, 
                            dr2_percent,
                            dr3_percent,
                            placebo_percent],
                           [dr1, dr2, dr3, placebo])
percent_change
```




    Capomulin   -0.194753
    Infubinol    0.461235
    Ketapril     0.570288
    Placebo      0.512980
    dtype: float64




```python
x_axis = np.arange(0, 4, 1)
barwidth = 0.8 

fig, ax = plt.subplots()

# masks for coloring bars
mask_minus = percent_change < 0
mask_plus = percent_change >= 0

rects1 = ax.bar(x_axis[mask_plus], percent_change[mask_plus], color="red", alpha=0.5, align="edge", width=barwidth)
rects2 = ax.bar(x_axis[mask_minus], percent_change[mask_minus], color="green", alpha=0.5, align="edge", width=barwidth)

# horizontal line at y=0, opacity = alpha
ax.hlines(0, -0.1, 10, alpha=0.5)
    
# limits of x-axis (with whitespace on left+right)
ax.set_xlim(-0.1, len(x_axis))

# limits of y-axis (with whitespace on top+bottom of bar chart)
ax.set_ylim(min(percent_change) - .1, max(percent_change) + .1)

# add some text for labels, title and axes ticks
ax.set_ylabel("% Tumor Volume Change", fontsize=14)
ax.set_xlabel("Treatment Drug", fontsize=14)
ax.set_title("% Tumor Change over 45-day Treatment Period", fontsize=18)

# ticks for bar chart's x-axis
ax.tick_params(direction="out", color="black", width=1, length=5, axis="y", pad=2)
ax.set_xticks([value+0.35 for value in x_axis])
ax.set_xticklabels(percent_change.index)

# padding below and to left of the ticks in the axes
ax.xaxis.labelpad = 25
ax.yaxis.labelpad = 25
ax.title.set_position([.5, 1.04])

# legend for the red vs. green coloring
ax.legend((rects1[0], rects2[0]), ('Tumor Growth', 'Tumor Reduction'), loc="upper left")

# new function!
def autolabel(rects):
    for rect in rects:
        height = rect.get_height() # for exact height of bar
        height_formatted = "{0:.0f}%".format(rect.get_height() * 100)  # format as a percentage for display
        ax.text(rect.get_x() + rect.get_width()/2., 0.5*height, height_formatted, ha='center', va='bottom')

autolabel(rects1)
autolabel(rects2)

plt.show()
```


![png](Images/output_32_0.png)

