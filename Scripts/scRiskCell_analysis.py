import pandas as pd
import numpy as np
import random
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score 
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.pyplot import MultipleLocator
from sklearn.preprocessing import label_binarize
from sklearn.metrics import roc_curve, auc
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

def PerformPCA(data,n_components):
    pca = PCA(n_components=n_components,random_state=42)
    pca_result = pca.fit_transform(data)
    pca_df = pd.DataFrame(pca_result, columns=[f'PC{i+1}' for i in range(n_components)])
    return pca_df

def DiseaseToLabel(disease):
    if disease == 'ND':
        return 0
    elif disease == 'preT2D':
        return 1
    elif disease == 'T2D':
        return 2

def PerformRegression(X,y):
    model = LogisticRegression(random_state=42)
    model.fit(X, y)
    return model

def GetDiaseaseIndex(state,model,pc_df_1,n_components,pc_df_2):
    if state == 0:
        index_1 = model.decision_function(pc_df_1.iloc[:,:n_components])
        index_2 = model.decision_function(pc_df_2.iloc[:,:n_components])
        index_df = pd.DataFrame({'disease index':np.concatenate((index_1, index_2)), 'disease':np.concatenate((pc_df_1.disease, pc_df_2.disease)),
                                 'donor':np.concatenate((pc_df_1.donor,pc_df_2.donor)),'label':np.concatenate((pc_df_1.label,pc_df_2.label)),
                                 'cell id':np.concatenate((pc_df_1.cell_id,pc_df_2.cell_id))})
        index_df = index_df.sort_values(by='label')
        return index_df
    elif state == 1:
        index = model.decision_function(pc_df_1.iloc[:,:n_components])
        index_df = pd.DataFrame({'disease index':index, 'disease':pc_df_1.disease,'donor':pc_df_1.donor,'label':pc_df_1.label,'cell id':pc_df_1.cell_id})
        index_df = index_df.sort_values(by='label')
        return index_df


def PlotVionlin(df,x,y,groups,colors,width,figsize):
    plt.figure(figsize=figsize)
    sns.violinplot(x=x, y=y, data=df,hue=groups, palette=colors,width=width)
    plt.xlabel('Disease')
    plt.ylabel('Disease Index')
    plt.gca().spines['top'].set_linewidth(2)
    plt.gca().spines['bottom'].set_linewidth(2)
    plt.gca().spines['left'].set_linewidth(2)
    plt.gca().spines['right'].set_linewidth(2)

def PlotBoxplot(df,x,y,groups,colors,width,figsize):
    plt.figure(figsize=figsize)
    sns.boxplot(x=x, y=y, data=df,hue=groups, palette=colors,width=width) 
    plt.xlabel('Disease')
    plt.ylabel('Disease Index')
    plt.gca().spines['top'].set_linewidth(2)
    plt.gca().spines['bottom'].set_linewidth(2)
    plt.gca().spines['left'].set_linewidth(2)
    plt.gca().spines['right'].set_linewidth(2)

def PlotDonorViolin(df,x,y,groups,colors,width,figsize):
    df.sort_values(by='label',ascending=True)
    plt.figure(figsize=figsize)
    sns.violinplot(x=x, y=y, data=df,palette=colors,width=width,inner="quart",hue=groups)
    plt.xticks(rotation=45, ha='right',va='top')
    plt.tight_layout()
    plt.legend(loc='upper left')
    plt.xlabel('Donor')
    plt.ylabel('Disease Index')
    plt.gca().spines['top'].set_linewidth(2)
    plt.gca().spines['bottom'].set_linewidth(2)
    plt.gca().spines['left'].set_linewidth(2)
    plt.gca().spines['right'].set_linewidth(2)

def PlotDonorBoxplot(df,x,y,groups,colors,width,figsize):
    df.sort_values(by='label',ascending=True)
    plt.figure(figsize=figsize)
    sns.boxplot(x=x, y=y,data=df,palette=colors,width=width,hue=groups)
    plt.xticks(rotation=45, ha='right',va='top')
    plt.tight_layout()
    plt.legend(loc='upper left')
    plt.xlabel('Donor')
    plt.ylabel('Disease Index')
    plt.gca().spines['top'].set_linewidth(2)
    plt.gca().spines['bottom'].set_linewidth(2)
    plt.gca().spines['left'].set_linewidth(2)
    plt.gca().spines['right'].set_linewidth(2)


# threshold 1 ：sort index quantile
def GetIndexValue_Quantile(df,ratio=0.85):
    df = df.sort_values(by='disease index',ascending=False)
    value = df['disease index'].quantile(ratio)
    return value

# threshold 2 ：slide windows
def GetIndexValue_SlideWindows(df,window_size=500,threshold=240,step=1):
    df = df.sort_values(by='disease index',ascending=False)
    count = 0
    start = 0
    value = None
    while start + window_size <= len(df):
        window = df['label'][start:start + window_size]
        count_2s = (window == 2).sum()
        if count_2s < threshold:
            break
        value = df['disease index'].iloc[start]
        start += step
        count += 1
    return count,start,value

def GetCellRisk(df,value):
    df['risk'] = (df['disease index'] > value).astype(int)
    return df

def GetDonorRiskRatio(df,data):
    risk_counts = df[df['risk'] == 1].groupby('donor').size()
    total_counts = df.groupby('donor').size()
    risk_ratios = (risk_counts / total_counts).fillna(0)
    risk_ratios = risk_ratios.reset_index()
    risk_ratios.columns = ['donor', 'risk ratios']
    risk_ratios['non risk ratios'] = 1 - risk_ratios['risk ratios']
    unique_data = data.drop_duplicates(subset=['donor'])
    risk_ratios['disease'] = risk_ratios['donor'].map(unique_data.set_index('donor')['disease'])
    risk_ratios['label'] = risk_ratios['disease'].apply(DiseaseToLabel)
    risk_ratios = risk_ratios.sort_values(by=['label','risk ratios'],ascending=[True,False])
    return risk_ratios

def GetSignificanceSymbol(p_value):
    if p_value < 0.001:
        return '***'
    elif p_value < 0.01:
        return '**'
    elif p_value < 0.05:
        return '*'
    else:
        return 'ns'
def AddStatAnnotation(ax, x1, x2, y, h, text):
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.8, c='k')
    plt.text((x1+x2)*.5, y+h, text, ha='center', va='bottom', color='k')


def PlotRatioBoxplot(state,df,x,y,groups,colors,width,linewidth,figsize,pYLoc1=0.8,pYLoc2=0.9,pYLoc3=1.0):
    plt.figure(figsize=figsize)
    plt.ylim(-0.05, 1.1)
    ax = sns.boxplot(x=x, y=y, data=df,hue=groups, palette=colors,width=width,linewidth=linewidth)
    if state == 0:
        group0 = df[df['label'].isin([0])]['risk ratios']
        group1 = df[df['label'].isin([1])]['risk ratios']
        group2 = df[df['label'].isin([2])]['risk ratios']
        _, p01 = mannwhitneyu(group0, group1, alternative="two-sided")
        _, p02 = mannwhitneyu(group0, group2, alternative="two-sided")
        _, p12 = mannwhitneyu(group1, group2, alternative="two-sided")
        
        _, p_adj, _, _ = multipletests([p01, p02, p12], method="fdr_bh")
        
        s01 = GetSignificanceSymbol(p_adj[0])
        s02 = GetSignificanceSymbol(p_adj[1])
        s12 = GetSignificanceSymbol(p_adj[2])
        
        AddStatAnnotation(ax, 0, 1, pYLoc1, 0.02, s01)
        AddStatAnnotation(ax, 0, 2, pYLoc3, 0.02, s02)
        AddStatAnnotation(ax, 1, 2, pYLoc2, 0.02, s12)

    if state == 1:
        group0 = df[df['label'].isin([0])]['risk ratios']
        group1 = df[df['label'].isin([2])]['risk ratios']
        _, p_val_0_1 = mannwhitneyu(group0, group1, alternative="two-sided")
        symbol_0_1 = GetSignificanceSymbol(p_val_0_1)
        AddStatAnnotation(ax, 0, 1, 1.0, 0.02, symbol_0_1)

    plt.xlabel('Disease')
    plt.ylabel('Risk cell ratios')
    plt.gca().spines['top'].set_linewidth(2)
    plt.gca().spines['bottom'].set_linewidth(2)
    plt.gca().spines['left'].set_linewidth(2)
    plt.gca().spines['right'].set_linewidth(2)


def PlotRatioStackBar(df,x,y,groups,colors,figsize):
    fig, ax = plt.subplots(figsize=figsize)
    bottom = np.zeros(len(df))
    for col, color in zip(['non risk ratios',  'risk ratios'], colors):
        ax.bar(df['donor'], df[col], bottom=bottom, label=col, color=color)
        bottom += df[col]
    ax.legend()
    ax.set_ylabel('Proportion')
    ax.set_xlabel('Donor')
    plt.xticks(rotation=45, ha='right',va='top')
    ax.margins(x=0.005)
    plt.tight_layout()
    plt.ylim(0.00, 1.02)
    plt.legend(loc='lower left')
    plt.gca().spines['top'].set_linewidth(2)
    plt.gca().spines['bottom'].set_linewidth(2)
    plt.gca().spines['left'].set_linewidth(2)
    plt.gca().spines['right'].set_linewidth(2)

def PlotBinaryROC(state,df,comparison_names,colors,figsize):
    if state == 0:
        df_01 = df[df['label'].isin([0, 1])]
        df_02 = df[df['label'].isin([0, 2])]
        df_12 = df[df['label'].isin([1, 2])]
        dfs = [df_01,df_02,df_12]
        dfs[1]['label'] = dfs[1]['label'].replace(2, 1)
        dfs[2]['label'] = dfs[2]['label'].replace(1, 0)
        dfs[2]['label'] = dfs[2]['label'].replace(2, 1)
        
        fpr = dict()
        tpr = dict()
        roc_auc = dict()
        
        for i,data in zip(range(3), dfs):
            fpr[i], tpr[i], _ = roc_curve(data['label'], data['risk ratios'])
            roc_auc[i] = auc(fpr[i], tpr[i])
        plt.figure(figsize=figsize)
        for i, color, name in zip(range(3), colors,comparison_names):
            plt.plot(fpr[i], tpr[i], color=color, lw=2,label='ROC curve of {0} (AUC = {1:0.3f})'.format(name, roc_auc[i]))
        plt.plot([0, 1], [0, 1], '--', lw=2,color='gray')
        plt.xlim([-0.05, 1.05])
        plt.ylim([-0.05, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('ROC curve ')
        plt.legend(loc="lower right")
        plt.grid(True,linestyle='--', linewidth=0.5, alpha=0.5)
        plt.gca().spines['top'].set_linewidth(2)
        plt.gca().spines['bottom'].set_linewidth(2)
        plt.gca().spines['left'].set_linewidth(2)
        plt.gca().spines['right'].set_linewidth(2)
    
    elif state == 1:
        df['label'] = df['label'].replace(2, 1)
        fpr, tpr, _ = roc_curve(df['label'], df['risk ratios'])
        roc_auc = auc(fpr, tpr)
        
        plt.figure(figsize=figsize)
        plt.plot(fpr, tpr, color=colors[0], lw=2,label='ROC curve (AUC = {0:0.3f})'.format(roc_auc))
        plt.plot([0, 1], [0, 1], '--', lw=2,color='gray')
        plt.xlim([-0.05, 1.05])
        plt.ylim([-0.05, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('ROC curve ')
        plt.legend(loc="lower right")
        plt.grid(True,linestyle='--', linewidth=0.5, alpha=0.5)
        plt.gca().spines['top'].set_linewidth(2)
        plt.gca().spines['bottom'].set_linewidth(2)
        plt.gca().spines['left'].set_linewidth(2)
        plt.gca().spines['right'].set_linewidth(2)    


def GetSVMClassifier(state,df):
    if state == 0:
        X, y = df[['risk ratios']], df.label
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
        svm_classifier = SVC(kernel='rbf', random_state=42,probability=True,decision_function_shape='ovr')
        svm_classifier.fit(X_train, y_train)
        return svm_classifier
    elif state == 1:
        X, y = df[['risk ratios']], df.label
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
        svm_classifier = SVC(kernel='rbf', random_state=42,probability=True)
        svm_classifier.fit(X_train, y_train)
        return svm_classifier

def GetLRClassifier(state,df):
    if state == 0:
        X, y = df[['risk ratios']], df.label
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
        lr_classifier = LogisticRegression(random_state=42,multi_class='multinomial')
        lr_classifier.fit(X_train, y_train)
        return lr_classifier
    elif state == 1:
        X, y = df[['risk ratios']], df.label
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
        lr_classifier = LogisticRegression(random_state=42)
        lr_classifier.fit(X_train, y_train)
        return lr_classifier

def PlotMultiROC(state,model,df,colors,title,figsize):
    if state == 0:
        X, y = df[['risk ratios']], df.label
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
        y_test_01 = label_binarize(y_test, classes=np.arange(3))
        fpr = {}
        tpr = {}
        roc_auc = {}
        for i in range(3):
            fpr[i], tpr[i], _ = roc_curve(y_test_01[:,i], model.predict_proba(X_test)[:, i])
            roc_auc[i] = auc(fpr[i], tpr[i])

        plt.figure(figsize=figsize)
        for i, color in zip(range(3), colors):
            plt.plot(fpr[i], tpr[i], color=color, lw=2,label=f'ROC curve of class {i} (AUC = {roc_auc[i]:0.3f})')
        plt.plot([0, 1], [0, 1], '--', lw=2,color='gray')
        plt.xlim([-0.05, 1.05])
        plt.ylim([-0.05, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title(title)
        plt.legend(loc="lower right")
        plt.grid(True,linestyle='--', linewidth=0.5, alpha=0.5)
        plt.gca().spines['top'].set_linewidth(2)
        plt.gca().spines['bottom'].set_linewidth(2)
        plt.gca().spines['left'].set_linewidth(2)
        plt.gca().spines['right'].set_linewidth(2)
    elif state == 1:
        X, y = df[['risk ratios']], df.label
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
        fpr, tpr, _ = roc_curve(y_test, model.predict_proba(X_test)[:, 1])
        roc_auc = auc(fpr, tpr)

        plt.figure(figsize=figsize)
        plt.plot(fpr, tpr, color=colors[0], lw=2,label=f'ROC curve (AUC = {roc_auc:0.3f})')
        plt.plot([0, 1], [0, 1], '--', lw=2,color='gray')
        plt.xlim([-0.05, 1.05])
        plt.ylim([-0.05, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title(title)
        plt.legend(loc="lower right")
        plt.grid(True,linestyle='--', linewidth=0.5, alpha=0.5)
        plt.gca().spines['top'].set_linewidth(2)
        plt.gca().spines['bottom'].set_linewidth(2)
        plt.gca().spines['left'].set_linewidth(2)
        plt.gca().spines['right'].set_linewidth(2)
