#!/usr/bin/env python
# coding: utf-8

# # Pymaceuticals Inc.
# ---
# 
# ### Analysis
# 
# - Add your analysis here.
#  

# In[232]:


# Dependencies and Setup
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st
from scipy.stats import linregress

# Study data files
mouse_metadata_path = "data/Mouse_metadata.csv"
study_results_path = "data/Study_results.csv"

# Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)

# Combine the data into a single dataset
data_df = pd.merge(study_results, mouse_metadata, on=["Mouse ID"], how = 'left')

# Display the data table for preview
data_df.head()


# In[183]:


# Checking the number of mice.
len(data_df["Mouse ID"].unique())


# In[184]:


# Getting the duplicate mice by ID number that shows up for Mouse ID and Timepoint. 
duplicates = data_df.loc[data_df.duplicated(subset = ['Mouse ID', 'Timepoint']), 'Mouse ID'].unique()
duplicates


# In[185]:


# Optional: Get all the data for the duplicate mouse ID. 
duplicates_list = data_df.loc[data_df['Mouse ID'] == 'g989']
duplicates_list


# In[186]:


# Create a clean DataFrame by dropping the duplicate mouse by its ID.
data_clean_df = data_df[~data_df['Mouse ID'].isin(duplicates)]
data_clean_df.head()


# In[187]:


# Checking the number of mice in the clean DataFrame.
len(data_clean_df["Mouse ID"].unique())


# ## Summary Statistics

# In[188]:


# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen
# Use groupby and summary statistical methods to calculate the following properties of each drug regimen: 
# mean, median, variance, standard deviation, and SEM of the tumor volume. 
# Assemble the resulting series into a single summary DataFrame.

mean = data_clean_df['Tumor Volume (mm3)'].groupby(data_clean_df['Drug Regimen']).mean()
median = data_clean_df['Tumor Volume (mm3)'].groupby(data_clean_df['Drug Regimen']).median()
variance = data_clean_df['Tumor Volume (mm3)'].groupby(data_clean_df['Drug Regimen']).var()
s_dev = data_clean_df['Tumor Volume (mm3)'].groupby(data_clean_df['Drug Regimen']).std()
sem = data_clean_df['Tumor Volume (mm3)'].groupby(data_clean_df['Drug Regimen']).sem()

summary_df = pd.DataFrame({"Mean Tumor Volume":mean, 
                           "Median Tumor Volume":median, 
                           "Tumor Volume Variance":variance, 
                           "Tumor Volume Std. Dev.":s_dev, 
                           "Tumor Volume Std. Err.":sem})
summary_df


# In[189]:


# Generate a summary statistics table of mean, median, variance, standard deviation, 
# and SEM of the tumor volume for each regimen

# Using the aggregation method, produce the same summary statistics in a single line.
summary_aggregation =  data_clean_df.groupby(['Drug Regimen'])[['Tumor Volume (mm3)']].agg(['mean', 'median', 'var', 'std', 'sem'])
summary_aggregation


# ## Bar and Pie Charts

# In[190]:


# Generate a bar plot showing the total number of timepoints for all mice tested for each drug regimen using Pandas.
plot_pd = data_clean_df["Drug Regimen"].value_counts().plot.bar()  
plt.xlabel("Drug Regimen")
plt.ylabel("Number of Mice Tested")
plt.show()


# In[191]:


# Generate a bar plot showing the total number of timepoints for all mice tested for each drug regimen using pyplot.
regimen_count = data_clean_df["Drug Regimen"].value_counts()
x_axis = regimen_count.index.values
y_axis = regimen_count.values

plt.bar(x_axis, y_axis, align='center')
plt.xlim(-0.75, len(x_axis)-0.25)
plt.xlabel("Drug Regimen")
plt.xticks(rotation="vertical")
plt.ylabel("Number of Mice Tested")

plt.show()


# In[192]:


# Generate a pie plot showing the distribution of female versus male mice using Pandas
gender_count = data_clean_df["Sex"].value_counts()
gender_count.plot.pie(autopct= "%1.1f%%")
plt.show()


# In[193]:


# Generate a pie plot showing the distribution of female versus male mice using pyplot
labels = ['Male', 'Female']
numbers = gender_count
plt.pie(numbers, labels=labels, autopct="%1.1f%%")
plt.ylabel('Sex')
plt.show()


# ## Quartiles, Outliers and Boxplots

# In[194]:


# Calculate the final tumor volume of each mouse across four of the treatment regimens:  
# Capomulin, Ramicane, Infubinol, and Ceftamin

# Start by getting the last (greatest) timepoint for each mouse
capomulin_time = data_clean_df.loc[data_clean_df["Drug Regimen"] == "Capomulin",:].groupby('Mouse ID').max()['Timepoint']
ramicane_time = data_clean_df.loc[data_clean_df["Drug Regimen"] == "Ramicane", :].groupby('Mouse ID').max()['Timepoint']
infubinol_time = data_clean_df.loc[data_clean_df["Drug Regimen"] == "Infubinol", :].groupby('Mouse ID').max()['Timepoint']
ceftamin_time = data_clean_df.loc[data_clean_df["Drug Regimen"] == "Ceftamin", :].groupby('Mouse ID').max()['Timepoint']

# Merge this group df with the original DataFrame to get the tumor volume at the last timepoint
capomulin_merged = pd.merge(capomulin_time, data_clean_df, on=("Mouse ID","Timepoint"),how="left")
ramicane_merged = pd.merge(ramicane_time, data_clean_df, on=("Mouse ID","Timepoint"),how="left")
infubinol_merged = pd.merge(infubinol_time, data_clean_df, on=("Mouse ID","Timepoint"),how="left")
ceftamin_merged = pd.merge(ceftamin_time, data_clean_df, on=("Mouse ID","Timepoint"),how="left")
infubinol_merged


# In[195]:


# Put treatments into a list for for loop (and later for plot labels)

# Create empty list to fill with tumor vol data (for plotting)
# Calculate the IQR and quantitatively determine if there are any potential outliers.   
    # Locate the rows which contain mice on each drug and get the tumor volumes  
    # add subset     
    # Determine outliers using upper and lower bounds
    
treatments = [capomulin_merged, ramicane_merged, infubinol_merged, ceftamin_merged]
tumor_vol_cap = capomulin_merged["Tumor Volume (mm3)"]
tumor_vol_ram = ramicane_merged["Tumor Volume (mm3)"]
tumor_vol_inf = infubinol_merged["Tumor Volume (mm3)"]
tumor_vol_cef = ceftamin_merged["Tumor Volume (mm3)"]
for treatment in treatments:
    quartiles = treatment["Tumor Volume (mm3)"].quantile([0.25, 0.5, 0.75])
    lowerq = quartiles[0.25]
    upperq = quartiles[0.75]
    iqr = upperq - lowerq
    lower_bound = lowerq - (iqr * 1.5)
    upper_bound = upperq + (iqr * 1.5)
    outlier = treatment[((treatment["Tumor Volume (mm3)"] < lower_bound) | (treatment["Tumor Volume (mm3)"] > upper_bound))]                  
    
    
    print(f"{treatment['Drug Regimen'][0]}'s potential outliers: {outlier['Tumor Volume (mm3)']}")
    #my number is not index 31 because I took them seperately


# In[196]:


# Generate a box plot that shows the distrubution of the tumor volume for each treatment group.
graph_data = [tumor_vol_cap, tumor_vol_ram, tumor_vol_inf, tumor_vol_cef]
regimen = ["Capomulin", "Ramicane", "Infubinol", "Ceftamin"]
plt.ylabel('Tumor Volume (mm3)')
plt.boxplot(graph_data, labels=regimen, vert=True)
plt.show()


# ## Line and Scatter Plots

# In[205]:


# Generate a line plot of tumor volume vs. time point for a mouse treated with Capomulin
specific_line = data_clean_df.loc[data_clean_df["Mouse ID"] == "l509"]
x_axis = specific_line["Timepoint"]
y_axis = specific_line["Tumor Volume (mm3)"]

plt.title('Capomulin treatmeant of mouse l509')
plt.plot(x_axis, y_axis)
plt.xlabel('Timepoint (days)')
plt.ylabel('Tumor Volume (mm3)')
plt.show()


# In[229]:


# Generate a scatter plot of average tumor volume vs. mouse weight for the Capomulin regimen
regi_cap = data_clean_df.loc[data_clean_df["Drug Regimen"] == "Capomulin"]
avg_regi_cap = regi_cap.groupby(["Mouse ID"]).mean()
x_axis = avg_regi_cap["Weight (g)"]
y_axis = avg_regi_cap["Tumor Volume (mm3)"]

plt.scatter(x_axis, y_axis)
plt.xlabel('Weight (g)')
plt.ylabel('Average Tumor Volume (mm3)')
plt.show()


# ## Correlation and Regression

# In[234]:


# Calculate the correlation coefficient and linear regression model 
# for mouse weight and average tumor volume for the Capomulin regimen
correlation = st.pearsonr(avg_regi_cap['Weight (g)'], avg_regi_cap['Tumor Volume (mm3)'])
print(f"The correlation between mouse weight and the average tumor volume is {round(correlation[0],2)}")

x_axis = avg_regi_cap["Weight (g)"]
y_axis = avg_regi_cap["Tumor Volume (mm3)"]
(slope, inter, rvalue, pvalue, stderr)= linregress(x_axis, y_axis)
fit = x_axis * slope + inter


plt.scatter(x_axis, y_axis)
plt.plot(x_axis, fit, color='red')
plt.xlabel("Weight (g)")
plt.ylabel("Average Tumor Volume (mm3)")
plt.show()


# In[ ]:




