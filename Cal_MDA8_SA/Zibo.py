'''
 # @ Author: Xue Jin
 # @ Create Time: 2021-07-05 18:52:31
 # @ Modified by: Xue Jin
 # @ Modified time: 2021-07-14 09:29:16
 # @ Description: 淄博项目专用函数库2.0
 '''


import pandas as pd
import numpy as np
import math
import re
from thefuzz import fuzz
from thefuzz import process
from chemicals import *
import json
from urllib import parse, request
import sys

def translate(sentence,src_lan='en',tgt_lan='zh',apikey='e6bc61094123576c9daf961192ffc61c'):
	url = 'https://api.niutrans.com/NiuTransServer/translation?'
	data = {"from":src_lan, "to":tgt_lan, "apikey":apikey, "src_text": sentence}
	data_en= parse.urlencode(data)
	req = url +"&"+ data_en
	res = request.urlopen(req)
	res = res.read()
	res_dict = json.loads(res)
	result=""
	if 'tgt_text' in res_dict:
		result = res_dict['tgt_text']
	else:
		result = res
	return result

def bar_labels(ax,offset=0.2,fontsize=12,fmt='%.1f',where='center',color='black'):
    if where=='up':
        for x in ax.containers:
            for y in x.get_children():
                ax.text(y.xy[0]+offset,y.get_height(),s=fmt%(y.get_height()),horizontalalignment='center',zorder=10,fontsize=fontsize,color=color)
        return ax
    elif where=='down':
        for x in ax.containers:
            for y in x.get_children():
                ax.text(y.xy[0]+offset,0,s=fmt%(y.get_height()),horizontalalignment='center',fontsize=fontsize,zorder=10,color=color)
        return ax
    elif where=='center':
        for x in ax.containers:
            for y in x.get_children():
                ax.text(y.xy[0]+offset,y.get_height(),s=fmt%(y.get_height()),horizontalalignment='center',va='top',fontsize=fontsize,zorder=10,color=color)
        return ax
    elif where=='top':
        ymax = ax.get_ylim()[1]
        for x in ax.containers:
            for y in x.get_children():
                ax.text(y.xy[0]+offset,ymax,s=fmt%(y.get_height()),horizontalalignment='center',va='bottom',fontsize=fontsize,zorder=10,color=color)
        return ax
    elif where=='bottom':
        ymin = ax.get_ylim()[0]
        for x in ax.containers:
            for y in x.get_children():
                ax.text(y.xy[0]+offset,ymin,s=fmt%(y.get_height()),horizontalalignment='center',va='top',fontsize=fontsize,zorder=10,color=color)
        return ax

def barh_labels(ax,offset=0.05,fontsize=12,fmt='%.1f',where='center',color='black'):
    if where=='top':
        for x in ax.containers:
            for y in x.get_children():
                ax.text(y.get_width(),y.xy[1]+offset,s=fmt%(y.get_width()),verticalalignment='center',zorder=10,fontsize=fontsize,color=color)
        return ax
    elif where=='bottom':
        for x in ax.containers:
            for y in x.get_children():
                ax.text(0,y.xy[1]+offset,s=fmt%(y.get_width()),verticalalignment='center',fontsize=fontsize,zorder=10,color=color)
        return ax
    elif where=='center':
        for x in ax.containers:
            for y in x.get_children():
                ax.text(y.get_width(),y.xy[1]+offset,s=fmt%(y.get_width()),verticalalignment='center',ha='right',fontsize=fontsize,zorder=10,color=color)
        return ax

def trans_season(x):
    if x.month in [12,1,2]:
        return '冬'
    elif x.month in [3,4,5]:
        return '春'
    elif x.month in [6,7,8]:
        return '夏'
    elif x.month in [9,10,11]:
        return '秋'

def rename_columns(data,newnames='CAS'):
    spec_info = pd.read_csv(r'D:\Desktop\codes\python codes\淄博\basefile\镇办站物种转化信息库_OBM.csv')
    spec_info = spec_info.set_index('species_name').to_dict('index')
    data.columns = [spec_info[columns][newnames] for columns in data.columns]
    return data

def name2cas(name):
    try:
        cas = CAS_from_any(name)
    except ValueError:
        cas = name
    return cas

def drop_othr(df):
    df.columns = [column.split('(')[0] for column in list(df.columns)]
    columns = list(df.columns)
    columns = [column.split('(')[0] for column in columns]
    df.columns = columns
    df = df.drop(columns=[
        '甲烷',
        '总烃',
        '非甲烷总烃',
        '甲硫醇',
        '乙硫醇',
        '二甲基硫',
        '二乙基硫',
        '二甲基二硫',
        '二硫化碳',
        '异丁烯',
        '1,3-丁二烯',
        '1,1,1-三氯乙烷',
        '1,1,2-三氯乙烷',
        '1,1-二氯乙烷',
        '1,1-二氯乙烯',
        '1,2-二氯乙烷',
        '1,2-二溴乙烷',
        '丙烯腈',
        '二氯甲烷',
        '氟里昂-114',
        '氟利昂12',
        '氯苯',
        '氯甲烷',
        '氯乙烷',
        '氯乙烯',
        '三氯甲烷',
        '顺-1,2-二氯乙烯',
        '四氯化碳',
        '四氯乙烯',
        '溴甲烷',
        '三氯乙烯',
        '四氯乙烯',
        '溴甲烷',
        '三氯乙烯',
        ])
    return df    

def cal_group(df_in,filter='no'):
    df_out = df_in.copy()
    df_out.columns = [column.split('(')[0] for column in list(df_out.columns)]
    VOClist = list(df_out.columns)
    alkane = [x for x in VOClist if ('烷' in x) and ('氯' not in x) and ('氟' not in x) and ('溴' not in x) and ('碘' not in x)]
    alkene = [x for x in VOClist if ('烯' in x) and ('氯' not in x) and ('氟' not in x) and ('溴' not in x) and ('碘' not in x) and ('苯' not in x)]  
    alkene.remove('异戊二烯')       
    bene = [x for x in VOClist if ('苯' in x) and ('氯' not in x) ]
    # df_out = df_out.reset_index()

    if filter=='yes':
        # df_out['烷烃'] = df_out[alkane].apply(lambda x: x.sum(min_count=1), axis=1)
        # df_out['烯烃'] = df_out[alkene + ['异戊二烯']].apply(lambda x: x.sum(min_count=1), axis=1)
        # df_out['芳香烃'] = df_out[bene].apply(lambda x: x.sum(min_count=1), axis=1)
        # df_out['VOCs'] = df_out[alkane + alkene + bene + ['乙炔']].apply(lambda x: x.sum(min_count=56//2), axis=1)
        # df_out = df_out[['烷烃','烯烃','芳香烃','乙炔','VOCs']]
        # cols = ['烷烃','烯烃','芳香烃','VOCs']
        # df_out[cols] = df_out[cols].replace({0:np.nan})

        df_out['烷烃'] = df_out[alkane].apply(lambda x: x.sum(min_count=1), axis=1)
        df_out['烯烃'] = df_out[alkene + ['异戊二烯']].apply(lambda x: x.sum(min_count=1), axis=1)
        df_out['芳香烃'] = df_out[bene].apply(lambda x: x.sum(min_count=1), axis=1)
        df_out['VOCs'] = df_out[alkane + alkene + bene].apply(lambda x: x.sum(min_count=56//2), axis=1)
        df_out = df_out[['烷烃','烯烃','芳香烃','VOCs']]
        cols = ['烷烃','烯烃','芳香烃','VOCs']
        df_out[cols] = df_out[cols].replace({0:np.nan})
        return df_out
    else:
        # df_out['烷烃'] = df_out[alkane].sum(axis=1,min_count=1)
        # df_out['烯烃'] = df_out[alkene].sum(axis=1,min_count=1)
        # df_out['芳香烃'] = df_out[bene].sum(axis=1,min_count=1)
        # df_out['VOCs'] = df_out[alkane + alkene + bene + ['异戊二烯'] + ['乙炔']].sum(axis=1,min_count=1)
        # df_out = df_out[['烷烃','烯烃','芳香烃','异戊二烯','乙炔','VOCs']]

        df_out['烷烃'] = df_out[alkane].sum(axis=1,min_count=1)
        df_out['烯烃'] = df_out[alkene].sum(axis=1,min_count=1)
        df_out['芳香烃'] = df_out[bene].sum(axis=1,min_count=1)
        df_out['VOCs'] = df_out[alkane + alkene + bene + ['异戊二烯']].sum(axis=1,min_count=1)
        df_out = df_out[['烷烃','烯烃','芳香烃','异戊二烯','VOCs']]
        return df_out

def get_group(df,content):
    df1 = df.copy()
    VOClist = list(df1.columns)
    alkane = [x for x in VOClist if ('烷' in x) and ('氯' not in x) and ('氟' not in x) and ('溴' not in x) and ('碘' not in x)]
    alkene = [x for x in VOClist if ('烯' in x) and ('氯' not in x) and ('氟' not in x) and ('溴' not in x) and ('碘' not in x) and ('苯' not in x)]  
    alkene.remove('异戊二烯')       
    bene = [x for x in VOClist if ('苯' in x) and ('氯' not in x) ]
    # df = df.reset_index()
    # df['烷烃'] = df[alkane].apply(lambda x: x.sum(), axis=1)
    # df['烯烃'] = df[alkene].apply(lambda x: x.sum(), axis=1)
    # df['芳香烃'] = df[bene].apply(lambda x: x.sum(), axis=1)
    # df['VOCs'] = df[alkane + alkene + bene + ['异戊二烯'] + ['乙炔']].apply(lambda x: x.sum(), axis=1)
    if content=='烯烃':
        df1 = df1[alkene]
    elif content=='芳香烃':
        df1 = df1[bene]
    return df1

def get_groupname(df):
    df1 = df.copy()
    VOClist = list(df1.columns)
    alkane = [x for x in VOClist if ('烷' in x) and ('氯' not in x) and ('氟' not in x) and ('溴' not in x) and ('碘' not in x)]
    alkene = [x for x in VOClist if ('烯' in x) and ('氯' not in x) and ('氟' not in x) and ('溴' not in x) and ('碘' not in x) and ('苯' not in x)]  
    alkene.remove('异戊二烯')       
    bene = [x for x in VOClist if ('苯' in x) and ('氯' not in x) ]
    # df = df.reset_index()
    # df['烷烃'] = df[alkane].apply(lambda x: x.sum(), axis=1)
    # df['烯烃'] = df[alkene].apply(lambda x: x.sum(), axis=1)
    # df['芳香烃'] = df[bene].apply(lambda x: x.sum(), axis=1)
    # df['VOCs'] = df[alkane + alkene + bene + ['异戊二烯'] + ['乙炔']].apply(lambda x: x.sum(), axis=1)
    return alkane,alkene,bene

def cal_ofp(df,type='cn'):
    df1 = df.copy()
    df1.columns = [column.split('(')[0] for column in list(df1.columns)]
    ofp = pd.read_csv(r'D:\Desktop\codes\python codes\淄博\basefile\镇办站物种转化信息库.csv',encoding='utf_8_sig')
    ofp = ofp.set_index('species_name')
    ofp = ofp.to_dict('index')
    VOC_spec = list(df1.columns)

    # for spec in VOC_spec:
    if type=='cn':
        for species in VOC_spec:
            if pd.isna(ofp[species]['MIR_cn']):
                df1[species] = df1[species]*ofp[species]['MIR_usa']
            else:
                df1[species] = df1[species]*ofp[species]['MIR_cn']
    elif type=='usa':
        for species in VOC_spec:
            df1[species] = df1[species]*ofp[species]['MIR_usa']
    return df1

def cal_vconc(df):
    df1 = df.copy()
    df1.columns = [column.split('(')[0] for column in list(df1.columns)]
    ofp = pd.read_csv(r'D:\Desktop\codes\python codes\淄博\basefile\镇办站物种转化信息库_OBM.csv',encoding='utf_8_sig')
    ofp = ofp.set_index('species_name')
    ofp = ofp.to_dict('index')
    VOC_spec = list(df1.columns)

    for spec in VOC_spec:
        # if (spec=='氮氧化物')|(spec=='NOx')|(spec=='PM2.5')|(spec=='气温')|(spec=='PM10')|(spec=='风速')|(spec=='风向')|(spec=='降水量')|(spec=='能见度')|(spec=='湿度')|(spec=='温度'):
        #     continue
        df1[spec] = df1[spec]/ofp[spec]['MWt']*22.4
    return df1

def cal_ppbc(df):
    df1 = df.copy()
    df1.columns = [column.split('(')[0] for column in list(df1.columns)]
    ofp = pd.read_csv(r'D:\Desktop\codes\python codes\淄博\basefile\镇办站物种转化信息库_OBM.csv',encoding='utf_8_sig')
    ofp = ofp.set_index('species_name')
    ofp['formula'] =  ofp['CAS'].apply(lambda x:search_chemical(x).formula)
    ofp['nC'] =  ofp['formula'].apply(lambda x:x[1])
    ofp['nC'] = ofp['nC'].astype(int)
    ofp = ofp.to_dict('index')
    VOC_spec = list(df1.columns)

    for spec in VOC_spec:
        if (spec=='氮氧化物')|(spec=='PM2.5')|(spec=='气温')|(spec=='PM10')|(spec=='风速')|(spec=='风向')|(spec=='降水量')|(spec=='能见度')|(spec=='湿度')|(spec=='温度'):
            continue
        df1[spec] = df1[spec]/ofp[spec]['MWt']*22.4*ofp[spec]['nC']
    return df1
    
def cal_soap(df):
    df1 = df.copy()
    df1.columns = [column.split('(')[0] for column in list(df1.columns)]
    ofp = pd.read_csv(r'D:\Desktop\codes\python codes\淄博\basefile\镇办站物种转化信息库.csv',encoding='utf_8_sig')
    ofp = ofp.set_index('species_name')
    ofp = ofp.to_dict('index')
    VOC_spec = list(df1.columns)

    for spec in VOC_spec:
        df1[spec] = ofp[spec]['soay']*df1[spec]
    return df1

def match_site(df,columnname,matchlist,paper=False):
    if paper:
        site_info = pd.read_csv(r'D:\Desktop\codes\python codes\淄博\basefile\淄博市镇办VOCs监测网络.csv',encoding='utf_8_sig')
    else:
        site_info = pd.read_csv(r'D:\Desktop\淄博\淄博市大气环境监测网络_10.21.csv',encoding='gbk')
    df = pd.merge(df,site_info[['站点名称'] + matchlist],how='left',left_on=columnname,right_on='站点名称')
    return df

    
def cal_ofpRank(df,num):
    df_groups = df.groupby(df.index)
    df_ofp = pd.DataFrame()
    for key,group in df_groups:
        group = group.T.sort_values(key,ascending=False)
        group = group.iloc[:num,:].reset_index()
        df_ofp = pd.concat([df_ofp,group],axis=1)
    return df_ofp


def trans_wd(x):
    wd = {'北':0,'北东北':22.5,'东北':45,'东东北':67.5,'东':90,'东东南':112.5,'东南':135,'南东南':157.5,'南':180,'南西南':202.5,'西南':225,'西西南':247.5,'西':270,'西西北':292.5,'西北':315,'北西北':337.5,np.nan:None}
    # if wd[x] True:
    x = wd[x]

    if x:
        y = (x)*math.pi/180
    else:
        y = np.nan

    return y

# def trans_wdAng(x):
#     y = x*math.pi/180
#     return y


def tranSpecName(tgtList):

    tgtList = {re.sub('[^\u4e00-\u9fa5^.^a-z^A-Z^0-9]','',spec):spec for spec in tgtList}

    species_info = pd.read_csv(r'D:\Desktop\codes\python codes\淄博\basefile\镇办站物种转化信息库.csv')
    species_info = species_info.set_index('species_name')
    eng2chn = species_info.to_dict('index')
    eng2chn = {value['species_Eng_name']:key for key,value in eng2chn.items()}

    PAMSname = species_info['species_Eng_name'].unique()
    PAMSname = {re.sub('[^\u4e00-\u9fa5^.^a-z^A-Z^0-9]','',spec):spec for spec in PAMSname}

    trans = {PAMSname[spec]:tgtList[process.extractBests(spec,tgtList.keys(),scorer=fuzz.token_sort_ratio,limit=1)[0][0]] \
            for spec in PAMSname.keys()}
    trans = pd.DataFrame(trans.values(),index=trans.keys())

    return trans

def boxplotFilterPro(col):
    # 计算iqr：数据四分之三分位值与四分之一分位值的差
    iqr = col.quantile(0.75)-col.quantile(0.25)
    # 根据iqr计算异常值判断阈值
    u_th = col.quantile(0.75) + 3*iqr # 上界
    l_th = col.quantile(0.25) - 3*iqr # 下界
    # 定义转换函数：如果数字大于上界则用上界值填充，小于下界则用下界值填充。
    def box_filter(x):
        if x > u_th:
            return np.nan
        elif x < l_th:
            return np.nan
        else:
            return x
    return col.map(box_filter)

def removeNegative(x):
    if x<0:
        return np.nan
    else:
        return x

def cleanData(VOCs,start,end):

    columnslist = VOCs.columns
    for column in columnslist:
        VOCs[column] = VOCs[column].apply(lambda x:removeNegative(x))

    #剔除大量零值行和大量缺值行
    VOCs = VOCs.drop(VOCs.loc[((VOCs == 0).astype(int).sum(axis=1) + VOCs.isna().sum(axis=1)>40)].index)

    #剔除重复行
    VOCs = VOCs.drop_duplicates(inplace=False)
    VOCs = VOCs.reset_index()

    timeindex = pd.date_range(start=start,end=end,freq='H')
    VOCs_groups = VOCs.groupby('站点名称')
    VOCs = pd.DataFrame()
    for key,group in VOCs_groups:
        group = group.set_index('时间')
        for column in columnslist:
            group[column] = boxplotFilterPro(group[column])
        group = group.reindex(timeindex)
        group['站点名称'] = key
        VOCs = VOCs.append(group)

    VOCs.index.name = '时间'
    VOCs = VOCs.reset_index()
    VOCs = VOCs.set_index(['时间','站点名称'])
    return VOCs

def OBMprep(data):
    data1 = cal_vconc(data)
    spec_info = pd.read_csv(r'D:\Desktop\codes\python codes\淄博\basefile\镇办站物种转化信息库.csv')
    spec_info = spec_info.set_index('species_name')
    spec_info = spec_info.to_dict('index')
    spec_info = {key:value['OBMname'] for key,value in spec_info.items()}
    data1 = data1.rename(columns=spec_info)
    data1['jno2'] = data1['SOLAR_R']*(8.87*(10**-6))+1.33*(10**-4)
    data1 = data1.reset_index()
    data1 = match_site(data1,'站点名称',['经度','纬度'])
    data1['YEAR'] = data1['时间'].dt.year
    data1['MONTH'] = data1['时间'].dt.month
    data1['DAY'] = data1['时间'].apply(lambda x:int(x.day))
    data1 = data1.rename(columns={'经度':'LON','纬度':'LAT','时间':'HOUR'})
    
    return data1