import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


#df_file = 'data/comparison_df.csv'
df_file = 'data/L1_comparison_df.csv'

df = pd.read_csv(df_file, index_col=0)
print(df)

#sns.lineplot(data=df, x='Size', y='Time', hue='Method')
sns.lineplot(data=df, x='Size', y='L1', hue='Method')
plt.show()