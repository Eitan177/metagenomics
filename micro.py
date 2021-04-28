import pdb
import numpy as np
import streamlit as st
from streamlit import caching 
#from bokeh.models.widgets import Button
#from bokeh.models import CustomJS
#from streamlit_bokeh_events import streamlit_bokeh_events
from io import StringIO
import pandas as pd
import altair as alt
import re
import SessionState


session = SessionState.get(run_id=0)
@st.cache(allow_output_mutation=True)
def get_df():
    return []

st.write('MetaGenomics Viewer')
new=st.button('Start a fresh session')
if new:
    caching.clear_cache()


df=get_df()

#copy_button = Button(label="Get Clipboard Data")

textuse=st.text_area('paste Clipboard')

#copy_button.js_on_event("button_click", CustomJS(code="""
#    navigator.clipboard.readText().then(text => document.dispatchEvent(new CustomEvent("GET_TEXT", {detail: text})))
#    """))
#result = streamlit_bokeh_events(
#    copy_button,
#    events="GET_TEXT",
#    key="get_text",
#    refresh_on_update=False,
#    override_height=75,
#    debounce_time=0)

if textuse:
    if True: #"GET_TEXT" in result:
        #dfnew = pd.read_table(StringIO(result.get("GET_TEXT")))
        dfnew = pd.read_table(StringIO(textuse))
        pasted=st.selectbox('Pasted from:',['choose which source','cDNA','DNA'],key=session.run_id)
        print(pasted)
        with st.beta_expander('show data'):
            st.table(dfnew)  
          
        if pasted !='choose which source':
            dfnew['source']=pasted
            df.append(dfnew)
            session.run_id += 1 


#pdb.set_trace()  


print(df)
#if len(df)>1 and np.all(df[-1]==df[-2]):
    #pdb.set_trace() and df[-2] != df[-1]:
#    session.run_id += 1

if  len(df)>0:
    for i in np.arange(0,len(df)):
        if df[i].shape[0]==0:
            del df[i]          
    
if len(df)==0:
    st.stop()


bound_df=pd.concat(df)

bound_df=bound_df.rename(columns={bound_df.columns[1]:'organism',bound_df.columns[5]:'reads'})


bound_df=bound_df.drop_duplicates()


def log2columnadd(df):
    
    df['reads']=df[['reads']].apply(lambda x: int(re.sub('[^0-9]','',str(x))),axis=1)
    df['readsLog2']=df[['reads']].apply(lambda x: np.log2(x),axis=1)
    return df   

df1=bound_df[['organism','reads','source']]
df1['kraken']='2016'
df2=bound_df.iloc[:,[9,13]]
df2['source']=df1['source']
df2['kraken']='2020'
df3=bound_df.iloc[:,[1,6]]
df3['source']=df1['source']
df3['kraken']='2016'

df4=bound_df.iloc[:,[9,14]]
df4['source']=df1['source']
df4['kraken']='2020'

df2.columns=df1.columns
df3.columns=df1.columns
df4.columns=df1.columns
df1=log2columnadd(df1)
df2=log2columnadd(df2)
df3=log2columnadd(df3)
df4=log2columnadd(df4)

controls=pd.concat((df1,df2),axis=0)
controls['sample']='NTC'

cases=pd.concat((df3,df4))
cases['sample']='sample'
casewcontrols=pd.concat((cases,controls),axis=0)
fofil=np.hstack((np.unique(casewcontrols['source']),np.unique(casewcontrols['kraken']),np.unique(casewcontrols['sample'])))

multifil=st.multiselect('remove',fofil,default=[])
distinct=st.selectbox('Distinguish',['sample type','kraken verison','source'])
distinctuse={'sample type':'sample','kraken verison':'kraken','source':'source'}

forp_distin=distinctuse[distinct]
casewcontrols=casewcontrols.reset_index()
sampfilt=[i  in multifil for i in casewcontrols['sample']]
sourcefilt=[i  in multifil for i in casewcontrols['source']]
krakenfilt=[i  in multifil for i in casewcontrols['kraken']]
remrows=np.where((pd.DataFrame((sampfilt,sourcefilt,krakenfilt))).T.apply(lambda x: any(x),axis=1))[0]

casewcontrols=casewcontrols.drop(index=remrows)
bb=alt.Chart(casewcontrols,width=1000).mark_bar().encode(x='readsLog2',y=alt.Y('organism', sort='-x'),color=forp_distin,tooltip=['organism','kraken','reads','source','sample']).interactive()
bb2=alt.Chart(casewcontrols).mark_bar().encode(x=alt.X('reads',scale=alt.Scale(domain=(0, int(controls[['reads']].quantile())))),y=alt.Y('organism', sort='-x'),color=forp_distin,tooltip=['organism','kraken','reads']).interactive()

with st.beta_expander('plots'):
    st.write(bb)


