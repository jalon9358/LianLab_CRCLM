import os
import scrublet as scr
import numpy as np
import pandas as pd

filename=os.listdir("./nonimmune/scrublet")

for i in range(0,len(filename)):
  os.chdir("./nonimmune/scrublet")
print(i)
data = pd.read_csv(filename[i],index_col = 0)
scrub = scr.Scrublet(data,expected_doublet_rate=0.05)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
#scrub.plot_histogram()
out=np.array([data.index,doublet_scores,predicted_doublets])
out_df=pd.DataFrame({'barcode':out[0,:],'score':out[1,:],'prediction':out[2,:]})
os.chdir("./nonimmune/doublet_score")
out_df.to_csv(f'{i+1}{".csv"}',index=False,header=True)
