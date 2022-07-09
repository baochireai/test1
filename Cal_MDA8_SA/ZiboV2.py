'''
 # @ Author: Xue Jin
 # @ Create Time: 2021-07-05 18:52:31
 # @ Modified by: Xue Jin
 # @ Modified time: 2021-08-18 17:45:16
 # @ Description: 淄博项目专用函数库2.0
 '''
import pandas as pd
import numpy as np
import math
import seaborn
import matplotlib.pyplot as plt

def vbar_label(ax,bar,fmt='%g'):
    for c in bar.containers:
        # Optional: if the segment is small or 0, customize the labels
        labels = [v.get_width() if v.get_width() > 0 else '' for v in c]
        # remove the labels parameter if it's not needed for customized labels
        ax.bar_label(c, labels=labels, label_type='center',fmt=fmt)

def trans_wd(x):
    wd = {'北':0,'北东北':22.5,'东北':45,'东东北':67.5,'东':90,'东东南':112.5,'东南':135,'南东南':157.5,'南':180,'南西南':202.5,'西南':225,'西西南':247.5,'西':270,'西西北':292.5,'西北':315,'北西北':337.5,np.nan:None}
    # if wd[x] True:
    x = wd[x]

    if x:
        y = (x)*math.pi/180
    else:
        y = np.nan

    return y


def match_site(df,columnname,matchlist):
    site_info = pd.read_csv(r'D:\Desktop\codes\python codes\淄博\basefile\淄博市大气环境监测网络_20220701.csv',encoding='utf_8_sig')
    df = pd.merge(df,site_info[['站点名称'] + matchlist],how='left',left_on=columnname,right_on='站点名称')
    return df

def trans_season(x):
    if x.month in [12,1,2]:
        return '冬'
    elif x.month in [3,4,5]:
        return '春'
    elif x.month in [6,7,8]:
        return '夏'
    elif x.month in [9,10,11]:
        return '秋'

class VOCzb:

    alkane = [
        '乙烷',
        '丙烷',
        '异丁烷',
        '正丁烷',
        '环戊烷',
        '异戊烷',
        '正戊烷',
        '2,2-二甲基丁烷',
        '2,3-二甲基丁烷',
        '2-甲基戊烷',
        '3-甲基戊烷',
        '正己烷',
        '甲基环戊烷',
        '2,4-二甲基戊烷',
        '环己烷',
        '2-甲基己烷',
        '2,3-二甲基戊烷',
        '3-甲基己烷',
        '2,2,4-三甲基戊烷',
        '正庚烷',
        '甲基环己烷',
        '2,3,4-三甲基戊烷',
        '2-甲基庚烷',
        '3-甲基庚烷',
        '正辛烷',
        '正壬烷',
        '正癸烷',
        '十一烷',
        '十二烷',]

    alkene = [
        '乙烯',
        '丙烯',
        '反-2-丁烯',
        '1-丁烯',
        '顺-2-丁烯',
        '反-2-戊烯',
        '1-戊烯',
        '顺-2-戊烯',
        '1-己烯',]

    arene = [
        '苯',
        '甲苯',
        '乙苯',
        '间,对二甲苯',
        '苯乙烯',
        '邻二甲苯',
        '异丙苯',
        '正丙苯',
        '间乙基甲苯',
        '对乙基甲苯',
        '1,3,5-三甲基苯',
        '1,2,4-三甲基苯',
        '1,2,3-三甲基苯',
        '邻乙基甲苯',
        '间二乙基苯',
        '对二乙基苯',]

    alkyne = ['乙炔']

    isop = ['异戊二烯']
    species = alkane + alkene + arene + alkyne + isop

    def __init__(self, path,time_col='时间',site_col='站点名称'):
        self.path = path
        self.data = pd.read_excel(path,parse_dates=[time_col])
        self.data.columns = [column.split('(')[0] for column in list(self.data.columns)]
        self.data = self.data[[time_col,site_col] + self.species]
        self.data = match_site(self.data,'站点名称',['区县'])
        self.data = self.data[[time_col,site_col,'区县'] + self.species]
        self.chem_info = pd.read_csv(r'D:\Desktop\codes\python codes\淄博\basefile\镇办站物种转化信息库.csv',encoding='utf_8_sig')
        self.chem_info = self.chem_info.set_index('species_name')
        self.chem_info = self.chem_info.to_dict('index')
    
    def get_histdata(self,start,end,sitenames=[]):
        self.data = pd.read_csv(r'D:\Desktop\上海大学工作成果汇总\淄博项目工作成果\2019-2021淄博市VOCs分析报告\VOCs案例库\2019-2021淄博市镇办VOCs数据_清洗后.csv',parse_dates=['时间'])
        self.data = self.data.set_index('时间')
        if len(sitenames) == 0:
            self.data = self.data
        elif len(sitenames) > 0:
            self.data = self.data.loc[self.data['站点名称'].isin([sitenames])]
        self.data = self.data.loc[start:end].reset_index()
        self.data = self.data[['时间','站点名称'] + self.species]
        self.data = match_site(self.data,'站点名称',['区县'])
        self.data = self.data[['时间','站点名称','区县'] + self.species]
        return self.data

    def cal_group(self,filter='no',index_list=['时间','区县','站点名称']):
        # df_out.columns = [column.split('(')[0] for column in list(df_out.columns)]

        if filter=='yes':
            self.data['烷烃'] = self.data[VOCzb.alkane].apply(lambda x: x.sum(min_count=1), axis=1)
            self.data['烯烃'] = self.data[VOCzb.alkene + ['异戊二烯']].apply(lambda x: x.sum(min_count=1), axis=1)
            self.data['芳香烃'] = self.data[VOCzb.arene].apply(lambda x: x.sum(min_count=1), axis=1)
            self.data['VOCs'] = self.data[self.arene + self.alkene + self.alkane].apply(lambda x: x.sum(min_count=56//2), axis=1)
            self.data = self.data[index_list + ['烷烃','烯烃','芳香烃','异戊二烯','乙炔','VOCs']]
            cols = ['烷烃','烯烃','芳香烃','VOCs']
            self.data[cols] = self.data[cols].replace({0:np.nan})
            return self.data
        else:
            self.data['烷烃'] = self.data[self.alkane].sum(axis=1,min_count=1)
            self.data['烯烃'] = self.data[self.alkene + ['异戊二烯']].sum(axis=1,min_count=1)
            self.data['芳香烃'] = self.data[self.arene].sum(axis=1,min_count=1)
            self.data['VOCs'] = self.data[self.arene + self.alkene + self.alkane].sum(axis=1,min_count=1)
            self.data = self.data[index_list + ['烷烃','烯烃','芳香烃','异戊二烯','乙炔','VOCs']]
            return self.data

    def cal_mean(self,by):
        self.data = self.data.groupby(by).mean()
        self.data = self.data.reset_index()
        return self.data

    def cal_ofp(self,type='usa'):
        if type=='cn':
            for species in self.species:
                if pd.isna(self.chem_info[species]['MIR_cn']):
                    self.data[species] = self.data[species]*self.chem_info[species]['MIR_usa']
                else:
                    self.data[species] = self.data[species]*self.chem_info[species]['MIR_cn']
        elif type=='usa':
            for species in self.species:
                self.data[species] = self.data[species]*self.chem_info[species]['MIR_usa']
        return self.data

    def cal_soap(self):
        for species in self.species:
            self.data[species] = self.data[species]*self.chem_info[species]['soay']
        return self.data

    def cal_vconc(self):
        for species in self.species:
            self.data[species] = self.data[species]/self.chem_info[species]['MWt']*22.4
        return self.data

    def rename_columns(self,newnames='OBMname'):
        spec_info = pd.read_csv(r'D:\Desktop\codes\python codes\淄博\basefile\镇办站物种转化信息库_OBM.csv')
        spec_info = spec_info.set_index('species_name').to_dict('index')
        spec_info = {key:value[newnames] for key,value in spec_info.items()}
        self.data.rename(columns=spec_info,inplace=True)
        return self.data

class VOCcz:
    
    # 加载超站VOCs物种信息库
    chem_info_path = r'D:\Desktop\codes\python codes\淄博\basefile\超站物种转化信息库.csv'
    chem_info = pd.read_csv(chem_info_path,encoding='utf_8_sig')
    pams_info = chem_info.loc[chem_info['type']=='PAMS']

    # 超站中各组分的名称
    species = list(chem_info['species_name'].values)
    alkanes = list(chem_info.loc[chem_info['Group']=='Alkanes']['species_name'].values)
    alkenes = list(chem_info.loc[chem_info['Group']=='Alkenes']['species_name'].values)
    aromatic_hydrocarbons = list(chem_info.loc[chem_info['Group']=='Aromatic_Hydrocarbons']['species_name'].values)
    ovoc = list(chem_info.loc[chem_info['Group']=='Oxygenated_Organics']['species_name'].values)
    halogenated_hydrocarbons = list(chem_info.loc[chem_info['Group']=='Halogenated_Hydrocarbons']['species_name'].values)

    # 超站中各组分的名称
    pams_species = list(pams_info['species_name'].values)
    pams_alkanes = list(pams_info.loc[pams_info['Group']=='Alkanes']['species_name'].values)
    pams_alkenes = list(pams_info.loc[pams_info['Group']=='Alkenes']['species_name'].values)
    pams_aromatic_hydrocarbons = list(pams_info.loc[pams_info['Group']=='Aromatic_Hydrocarbons']['species_name'].values)
    pams_ovoc = list(pams_info.loc[pams_info['Group']=='Oxygenated_Organics']['species_name'].values)
    pams_halogenated_hydrocarbons = list(pams_info.loc[pams_info['Group']=='Halogenated_Hydrocarbons']['species_name'].values)

    # 创建超站VOCs物种信息库索引
    chem_info_index = chem_info.copy()
    chem_info_index = chem_info_index.set_index('species_name')
    chem_info_index = chem_info_index.to_dict('index')


    def __init__(self, path,time_col='采集时间'):
        self.data = pd.read_html(path)[0]
        self.data[time_col] = pd.to_datetime(self.data[time_col])
        # self.data = pd.read_excel(path,parse_dates=[time_col])
        self.data = self.data[[time_col] + self.species]

    def cal_group(self,filter='no',index_list=['采集时间'],inplace='yes'):
        data = self.data.copy()
        if filter=='yes':
            data['烷烃'] = data[self.alkanes].apply(lambda x: x.sum(min_count=1), axis=1)
            data['烯烃'] = data[self.alkenes + ['异戊二烯']].apply(lambda x: x.sum(min_count=1), axis=1)
            data['芳香烃'] = data[self.aromatic_hydrocarbons].apply(lambda x: x.sum(min_count=1), axis=1)
            data['VOCs'] = data[self.alkanes + self.alkenes + self.aromatic_hydrocarbons + ['乙炔','异戊二烯']].apply(lambda x: x.sum(min_count=56//2), axis=1)
            data['卤代烃'] = data[self.halogenated_hydrocarbons].apply(lambda x: x.sum(min_count=1), axis=1)
            data['OVOC'] = data[self.ovoc].apply(lambda x: x.sum(min_count=1), axis=1)
            data = data[index_list + ['烷烃','烯烃','芳香烃','乙炔','异戊二烯','VOCs','卤代烃','OVOC']]
            cols = ['烷烃','烯烃','芳香烃','乙炔','VOCs','卤代烃','OVOC']
            data[cols] = data[cols].replace({0:np.nan})

        elif filter=='no':
            data['烷烃'] = data[self.alkanes].apply(lambda x: x.sum(min_count=1), axis=1)
            data['烯烃'] = data[self.alkenes + ['异戊二烯']].apply(lambda x: x.sum(min_count=1), axis=1)
            data['芳香烃'] = data[self.aromatic_hydrocarbons].apply(lambda x: x.sum(min_count=1), axis=1)
            data['VOCs'] = data[self.alkanes + self.alkenes + self.aromatic_hydrocarbons + ['乙炔','异戊二烯']].apply(lambda x: x.sum(min_count=56//2), axis=1)
            data['卤代烃'] = data[self.halogenated_hydrocarbons].apply(lambda x: x.sum(min_count=1), axis=1)
            data['OVOC'] = data[self.ovoc].apply(lambda x: x.sum(min_count=1), axis=1)
            data = data[index_list + ['烷烃','烯烃','芳香烃','乙炔','异戊二烯','VOCs','卤代烃','OVOC']]
        if inplace=='yes':
            self.data = data
            return self.data
        elif inplace=='no':
            return data

    def cal_pams_group(self,filter='no',index_list=['采集时间'],inplace='yes'):
        data = self.data.copy()
        if filter=='yes':
            data['烷烃'] = data[self.pams_alkanes].apply(lambda x: x.sum(min_count=1), axis=1)
            data['烯烃'] = data[self.pams_alkenes + ['异戊二烯']].apply(lambda x: x.sum(min_count=1), axis=1)
            data['芳香烃'] = data[self.pams_aromatic_hydrocarbons].apply(lambda x: x.sum(min_count=1), axis=1)
            data['VOCs'] = data[self.pams_alkanes + self.pams_alkenes + self.pams_aromatic_hydrocarbons + ['乙炔','异戊二烯']].apply(lambda x: x.sum(min_count=56//2), axis=1)
            data = data[index_list + ['烷烃','烯烃','芳香烃','乙炔','异戊二烯','VOCs']]
            cols = ['烷烃','烯烃','芳香烃','乙炔','VOCs']
            data[cols] = data[cols].replace({0:np.nan})

        elif filter=='no':
            data['烷烃'] = data[self.pams_alkanes].apply(lambda x: x.sum(min_count=1), axis=1)
            data['烯烃'] = data[self.pams_alkenes + ['异戊二烯']].apply(lambda x: x.sum(min_count=1), axis=1)
            data['芳香烃'] = data[self.pams_aromatic_hydrocarbons].apply(lambda x: x.sum(min_count=1), axis=1)
            data['VOCs'] = data[self.pams_alkanes + self.pams_alkenes + self.pams_aromatic_hydrocarbons + ['乙炔','异戊二烯']].apply(lambda x: x.sum(min_count=56//2), axis=1)
            data = data[index_list + ['烷烃','烯烃','芳香烃','乙炔','异戊二烯','VOCs']]
        if inplace=='yes':
            self.data = data
            return self.data
        elif inplace=='no':
            return data

    def cal_ofp(self,type='usa',unit='ppbv',inplace='yes'):
        data = self.data.copy()
        data = self.cal_mconc(inplace='no')
        if unit=='ppbv':
            if type=='cn':
                for species in self.species:
                    if pd.isna(self.chem_info_index[species]['MIR_cn']):
                        data[species] = data[species]*self.chem_info_index[species]['MIR_usa']/48*22.4
                    else:
                        data[species] = data[species]*self.chem_info_index[species]['MIR_cn']/48*22.4
            elif type=='usa':
                for species in self.species:
                    data[species] = data[species]*self.chem_info_index[species]['MIR_usa']/48*22.4
        elif unit=='ugm3':
            if type=='cn':
                for species in self.species:
                    if pd.isna(self.chem_info_index[species]['MIR_cn']):
                        data[species] = data[species]*self.chem_info_index[species]['MIR_usa']
                    else:
                        data[species] = data[species]*self.chem_info_index[species]['MIR_cn']
            elif type=='usa':
                for species in self.species:
                    data[species] = data[species]*self.chem_info_index[species]['MIR_usa']
        if inplace=='yes':
            self.data = data
            return self.data
        elif inplace=='no':
            return data

    def cal_mconc(self,inplace='yes'):
        data = self.data.copy()
        for species in self.species:
                data[species] = data[species]*self.chem_info_index[species]['MWt']/22.4
        if inplace=='yes':
            self.data = data
            return self.data
        elif inplace=='no':
            return data