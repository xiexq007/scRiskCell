import pandas as pd
import numpy as np
import random
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
import seaborn as sns

def Split(df, category_group1, category_group2):
    endpoint_df = df[df['Category'].isin(category_group1)].reset_index(drop=True)
    intermediate_df = df[df['Category'].isin(category_group2)].reset_index(drop=True)
    return endpoint_df, intermediate_df


def PerformPCA(endpoint_df, intermediate_df, n_components=20):
    data_for_pca_ep = endpoint_df.iloc[:, :-4]
    pca_ep = PCA(n_components=n_components, random_state=42)
    pca_result_ep = pca_ep.fit_transform(data_for_pca_ep)
    pca_df_ep = pd.DataFrame(pca_result_ep, columns=[f'PC{i+1}' for i in range(n_components)])
    last4_cols_ep = endpoint_df.iloc[:, -4:].reset_index(drop=True)
    pca_endpoint_df = pd.concat([pca_df_ep, last4_cols_ep], axis=1)

    if intermediate_df.empty:
        pca_intermediate_df = pd.DataFrame()
    else:
        data_for_pca_int = intermediate_df.iloc[:, :-4]
        pca_int = PCA(n_components=n_components, random_state=42)
        pca_result_int = pca_int.fit_transform(data_for_pca_int)
        pca_df_int = pd.DataFrame(pca_result_int, columns=[f'PC{i+1}' for i in range(n_components)])
        last4_cols_int = intermediate_df.iloc[:, -4:].reset_index(drop=True)
        pca_intermediate_df = pd.concat([pca_df_int, last4_cols_int], axis=1)

    return pca_endpoint_df, pca_intermediate_df


def PerformRegression(pca_endpoint_df, n_components=20):
    model = LogisticRegression(random_state=42)
    X = pca_endpoint_df.iloc[:, :n_components]    
    y = pca_endpoint_df['Label']         
    model.fit(X, y)
    return model


def GetDiaseaseIndex(pca_endpoint_df, pca_intermediate_df, model, n_components=20):
    if pca_intermediate_df.empty:
        index = model.decision_function(pca_endpoint_df.iloc[:,:n_components])
        index_df = pd.DataFrame({'Donor':pca_endpoint_df.Donor, 'Category':pca_endpoint_df.Category, 'Label':pca_endpoint_df.Label, 'Cell_id':pca_endpoint_df.Cell_id, 'Disease_index':index})
        index_df = index_df.sort_values(by='Label')
        return index_df
    else:
        index_ep = model.decision_function(pca_endpoint_df.iloc[:,:n_components])
        index_int = model.decision_function(pca_intermediate_df.iloc[:,:n_components])
        index_df = pd.DataFrame({'Donor':np.concatenate((pca_endpoint_df.Donor, pca_intermediate_df.Donor)), 
                                 'Category':np.concatenate((pca_endpoint_df.Category, pca_intermediate_df.Category)),
                                 'Label':np.concatenate((pca_endpoint_df.Label, pca_intermediate_df.Label)),
                                 'Cell_id':np.concatenate((pca_endpoint_df.Cell_id, pca_intermediate_df.Cell_id)),
                                 'Disease_index':np.concatenate((index_ep, index_int))})
        index_df = index_df.sort_values(by='Label')
        return index_df


def GetIndexThreshold(index_df, get_threshold_params):
    if not get_threshold_params:
        ratio = 0.85
        sorted_df = index_df.sort_values(by='Disease_index', ascending=False)
        value = sorted_df['Disease_index'].quantile(ratio)
        return value
    elif len(get_threshold_params) == 3:
        window_size = get_threshold_params.get('window_size')
        threshold = get_threshold_params.get('threshold')
        step = get_threshold_params.get('step')
        max_label = index_df['Label'].max()
        sorted_df = index_df.sort_values(by='Disease_index', ascending=False)
        start = 0
        value = None
        while start + window_size <= len(sorted_df):
            window = sorted_df['Label'][start:start + window_size]
            counts = (window == max_label).sum()
            if counts < threshold:
                break
            value = sorted_df['Disease_index'].iloc[start]
            start += step
        return value
    else:
        ratio = get_threshold_params.get('ratio')
        sorted_df = index_df.sort_values(by='Disease_index', ascending=False)
        value = sorted_df['Disease_index'].quantile(ratio)
        return value


def GetCellRisk(index_df, threshold):
    risk_df = index_df.copy()  
    risk_df['Risk'] = (risk_df['Disease_index'] > threshold).astype(int)
    return risk_df


def GetDonorRiskRatio(risk_df):
    risk_counts = risk_df[risk_df['Risk'] == 1].groupby('Donor').size()
    total_counts = risk_df.groupby('Donor').size()
    risk_ratios = (risk_counts / total_counts).fillna(0)
    risk_ratios = risk_ratios.reset_index()
    risk_ratios.columns = ['Donor', 'Risk_ratio']
    risk_ratios['non_Risk_ratio'] = 1 - risk_ratios['Risk_ratio']
    unique_data = risk_df.drop_duplicates(subset=['Donor'])
    risk_ratios['Category'] = risk_ratios['Donor'].map(unique_data.set_index('Donor')['Category'])
    risk_ratios['Label'] = risk_ratios['Donor'].map(unique_data.set_index('Donor')['Label'])
    risk_ratios = risk_ratios[['Donor', 'Category', 'Label', 'Risk_ratio', 'non_Risk_ratio']]
    risk_ratios = risk_ratios.sort_values(by=['Label','Risk_ratio'],ascending=[True, False])
    return risk_ratios


def scRiskCell(df, category_group1, category_group2, n_components=20, get_threshold_params=None):
    endpoint_df, intermediate_df = Split(df, category_group1, category_group2)
    pca_endpoint_df, pca_intermediate_df = PerformPCA(endpoint_df, intermediate_df, n_components)
    model = PerformRegression(pca_endpoint_df, n_components)
    index_df = GetDiaseaseIndex(pca_endpoint_df, pca_intermediate_df, model, n_components)
    threshold = GetIndexThreshold(index_df, get_threshold_params)
    risk_df = GetCellRisk(index_df, threshold)
    donor_risk_ratio = GetDonorRiskRatio(risk_df)

    return pca_endpoint_df, pca_intermediate_df, model, index_df, threshold, risk_df, donor_risk_ratio


def PlotIndexViolin_by_Category(index_df, colors=None, save_path=None):
    n_colors = index_df['Label'].nunique()
    if colors is None:
        colors = sns.color_palette('Set2', n_colors).as_hex()
    width = max(5, 3 + n_colors * 1)
    plt.figure(figsize=(width, 5))
    sns.violinplot(x='Category', y='Disease_index', data=index_df, hue='Category', palette=colors, width=0.65)
    plt.xlabel('Category')
    plt.ylabel('Disease index')
    plt.gca().spines['top'].set_linewidth(1.5)
    plt.gca().spines['bottom'].set_linewidth(1.5)
    plt.gca().spines['left'].set_linewidth(1.5)
    plt.gca().spines['right'].set_linewidth(1.5)

    if save_path is not None:
        plt.savefig(save_path, bbox_inches='tight')
    plt.show()


def PlotIndexBoxplot_by_Category(index_df, colors=None, save_path=None):
    n_colors = index_df['Label'].nunique()
    if colors is None:
        colors = sns.color_palette('Set2', n_colors).as_hex()
    width = max(5, 3 + n_colors * 1)
    plt.figure(figsize=(width, 5))
    sns.boxplot(x='Category', y='Disease_index', data=index_df, hue='Category', palette=colors, width=0.65)
    plt.xlabel('Category')
    plt.ylabel('Disease index')
    plt.gca().spines['top'].set_linewidth(1.5)
    plt.gca().spines['bottom'].set_linewidth(1.5)
    plt.gca().spines['left'].set_linewidth(1.5)
    plt.gca().spines['right'].set_linewidth(1.5)

    if save_path is not None:
        plt.savefig(save_path, bbox_inches='tight')
    plt.show()

def PlotIndexViolin_by_Donor(index_df, colors=None, save_path=None):
    n_colors = index_df['Label'].nunique()
    if colors is None:
        colors = sns.color_palette('Set2', n_colors).as_hex()
    n_donors = index_df['Donor'].nunique()
    width = max(5, 0.5 + n_donors * 0.3)
    plt.figure(figsize=(width, 5))
    sns.violinplot(x='Donor', y='Disease_index', data=index_df, hue='Category', palette=colors, width=0.85, inner="quart")
    plt.xticks(rotation=45, ha='right',va='top')
    plt.tight_layout()
    plt.legend(loc='upper left')
    plt.xlabel('Donor')
    plt.ylabel('Disease index')
    plt.gca().spines['top'].set_linewidth(1.5)
    plt.gca().spines['bottom'].set_linewidth(1.5)
    plt.gca().spines['left'].set_linewidth(1.5)
    plt.gca().spines['right'].set_linewidth(1.5)

    if save_path is not None:
        plt.savefig(save_path, bbox_inches='tight')
    plt.show()


def PlotIndexBoxplot_by_Donor(index_df, colors=None, save_path=None):
    n_colors = index_df['Label'].nunique()
    if colors is None:
        colors = sns.color_palette('Set2', n_colors).as_hex()
    n_donors = index_df['Donor'].nunique()
    width = max(5, 0.5 + n_donors * 0.3)
    plt.figure(figsize=(width, 5))
    sns.boxplot(x='Donor', y='Disease_index', data=index_df, hue='Category', palette=colors, width=0.85)
    plt.xticks(rotation=45, ha='right',va='top')
    plt.tight_layout()
    plt.legend(loc='upper left')
    plt.xlabel('Donor')
    plt.ylabel('Disease index')
    plt.gca().spines['top'].set_linewidth(1.5)
    plt.gca().spines['bottom'].set_linewidth(1.5)
    plt.gca().spines['left'].set_linewidth(1.5)
    plt.gca().spines['right'].set_linewidth(1.5)

    if save_path is not None:
        plt.savefig(save_path, bbox_inches='tight')
    plt.show()


def PlotRatioBoxplot_by_Category(donor_risk_ratio, colors=None, save_path=None):
    n_colors = donor_risk_ratio['Label'].nunique()
    if colors is None:
        colors = sns.color_palette('Set2', n_colors).as_hex()
    width = max(5, 3 + n_colors * 1)
    plt.figure(figsize=(width, 5))
    plt.ylim(-0.03, 1.03)
    sns.boxplot(x='Category', y='Risk_ratio', data=donor_risk_ratio, hue='Category', palette=colors, width=0.65)
    plt.xlabel('Category')
    plt.ylabel('Risk cell ratio')
    plt.gca().spines['top'].set_linewidth(1.5)
    plt.gca().spines['bottom'].set_linewidth(1.5)
    plt.gca().spines['left'].set_linewidth(1.5)
    plt.gca().spines['right'].set_linewidth(1.5)

    if save_path is not None:
        plt.savefig(save_path, bbox_inches='tight')
    plt.show()


def PlotRatioStackBar_by_Donor(donor_risk_ratio, colors=None, save_path=None):
    if colors is None:
        colors = sns.color_palette('Set2', 2).as_hex()
    n_donors = donor_risk_ratio['Donor'].nunique()
    width = max(1.5, 0.5 + n_donors * 0.3)
    figsize=(width, 5)
    fig, ax = plt.subplots(figsize=figsize)
    bottom = np.zeros(len(donor_risk_ratio))
    for col, color in zip(['non_Risk_ratio',  'Risk_ratio'], colors):
        ax.bar(donor_risk_ratio['Donor'], donor_risk_ratio[col], bottom=bottom, label=col, color=color)
        bottom += donor_risk_ratio[col]
    ax.legend()
    ax.set_ylabel('Proportion')
    ax.set_xlabel('Donor')
    plt.xticks(rotation=45, ha='right',va='top')
    ax.margins(x=0.005)
    plt.tight_layout()
    plt.ylim(0.00, 1.02)
    plt.legend(loc='lower left')
    plt.gca().spines['top'].set_linewidth(1.5)
    plt.gca().spines['bottom'].set_linewidth(1.5)
    plt.gca().spines['left'].set_linewidth(1.5)
    plt.gca().spines['right'].set_linewidth(1.5)
    
    if save_path is not None:
        plt.savefig(save_path, bbox_inches='tight')
    plt.show()


def PlotROC(donor_risk_ratio, color=None, save_path=None):
    max_label = donor_risk_ratio['Label'].max()
    new_df = donor_risk_ratio[donor_risk_ratio['Label'].isin([0, max_label])].copy()
    new_df['Label'] = new_df['Label'].replace(max_label, 1)
    fpr, tpr, _ = roc_curve(new_df['Label'], new_df['Risk_ratio'])
    roc_auc = auc(fpr, tpr)
    if color is None:
        color = sns.color_palette('Set2', 1).as_hex()[0]
    plt.figure(figsize=(5,5))
    plt.plot(fpr, tpr, color=color, lw=2, label='ROC curve (AUC = {0:0.3f})'.format(roc_auc))
    plt.plot([0, 1], [0, 1], '--', lw=2, color='gray')
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC curve ')
    plt.legend(loc="lower right")
    plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.5)
    plt.gca().spines['top'].set_linewidth(1.5)
    plt.gca().spines['bottom'].set_linewidth(1.5)
    plt.gca().spines['left'].set_linewidth(1.5)
    plt.gca().spines['right'].set_linewidth(1.5)    

    if save_path is not None:
        plt.savefig(save_path, bbox_inches='tight')
    plt.show()