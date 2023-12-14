#!/usr/bin/env python

import numpy as np
import pandas as pd
from pydeseq2 import preprocessing
from matplotlib import pyplot as plt



# Step 2.1: Initial exploration=========================================================================================================================

# Using Python (pandas, numpy, matplotlib, etc.) explore these data with your partner, searching for any interesting features or patterns. 
# Jot down any interesting patterns you observe as notes (no need to submit). For each feature/pattern you observe, think about what kind of plot would best communicate 
# that feature/pattern.


df=pd.read_csv("life_expectancy.csv", index_col=0)
print(df)

#life expectancy (LE) by entity
mean_LE_by_country=df.groupby('Entity')['LifeExpectancy'].mean()
#print(mean_LE_by_country)

#countries= mean["Code"]

if 'Code' in df.columns:
    entity_column = df['Code']
    # Your code here using 'entity_column'
else:
    print("Column 'Entity' does not exist.")




#LE by year

mean_LE_by_year=df.groupby('Year').agg({'LifeExpectancy' : 'mean'})
#print(mean_LE_by_year)




#LE by year and entity
#print(df)
df.reset_index(inplace=True)
df_grouped = df.groupby(['Entity','Year'])
#print(df_grouped.first())

mean_country_year=df.groupby(['Entity','Year'])['LifeExpectancy'].mean()
# mean_country_year_df=pd.DataFrame(mean_country_year)


#print(mean_country_year)






#Mean year of measurement for each entity
mean_year_by_country=df.groupby('Entity')['Year'].mean()
#print(mean_year_by_country)


#print(mean_year_by_country["Italy"])

#global mean LE over all time
mean_LE_all= df["LifeExpectancy"].mean()
#print(mean_LE_all)


# print(mean_LE_by_country)
# print(type(mean_LE_by_country


low_LE={}
for i in df.index:
    if df.loc[i, "LifeExpectancy"] < 30:

        low_LE.setdefault(df.loc[i, "Entity"], 0)
        low_LE[df.loc[i, "Entity"]] += 1

#print(low_LE)


high_LE={}
for i in df.index:
    if df.loc[i, "LifeExpectancy"] > 83:

        high_LE.setdefault(df.loc[i, "Entity"], 0)
        high_LE[df.loc[i, "Entity"]] += 1

print(high_LE)

#STEP2.2: plotting the trends =============================================================================
fig, ax=plt.subplots()

mean_LE_by_country = mean_LE_by_country.sort_values()

ax.scatter(mean_LE_by_country.index, mean_LE_by_country)
ax.set_xticklabels(labels = mean_LE_by_country.index, rotation=90, size = 3)
plt.title('Mean Life Expectancy for each Entity')
plt.xlabel('Entity')
plt.ylabel('Life Expectancy (Years)')


#mean_LE_by_year.plot("LifeExpectancy", kind="scatter")
#plt.title('Mean Life Expectancy per Year for each Country')
#plt.xlabel('Year')
#plt.ylabel('Values')
#plt.legend()

fig, ax=plt.subplots()
ax.scatter(mean_LE_by_year.index, mean_LE_by_year)
plt.title('Mean Global Life Expectancy per Year')
plt.xlabel('Year')
plt.ylabel('Life Expectancy (Years)')



fig, ax=plt.subplots()
ax.scatter(low_LE.keys(), low_LE.values())
ax.set_xticklabels(labels = low_LE.keys(), rotation=90, size = 7)
plt.title('Countries with life expectancy below 30')
plt.ylabel('Number of years with life expectancy below 30 years')
plt.xlabel('Entity')





plt.tight_layout()
plt.show()



