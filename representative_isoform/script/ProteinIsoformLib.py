import numpy as np


def average_PSM(row, df):
    header = list(df.columns)
    select_PSMs=[]
    for names in header:
        if "PSM" in names:
            select_PSMs.append(names)
    series = row[select_PSMs].astype('float')
#     print (select_PSMs)
    
    avg_psm = series.sum(skipna=True)/len(select_PSMs)
    return avg_psm



def average_intensity(row, df):
    header = list(df.columns)
    select_channel=[]
    for names in header:
        if "sig" in names:
            select_channel.append(names)
#     avg_intensity = row[select_channel].sum(skipna=True)/len(select_channel)
    series = row[select_channel].astype('float')
#     print (select_channel)
    avg_intensity = series.sum(skipna=True)/len(select_channel)
    return avg_intensity



def get_repeat_index(list1, max_value):
    index_list = []
    for index, batch_cnt in enumerate(list1):
        if batch_cnt == max_value:
            index_list.append(index)
    return index_list
            

def select_correct_iso_lists(access_list, batch_index_list, intensity_list_avg, psm_list_avg):
    access_list_new = []
    intensity_list_new = []
    psm_list_new = []
    
    for values in batch_index_list:
        access_list_new.append(access_list[values])
        intensity_list_new.append(intensity_list_avg[values])
        psm_list_new.append(psm_list_avg[values])
    
    return access_list_new, intensity_list_new, psm_list_new



'''This function looks for the priority sp| signature in the uniprot dataset and generate a unique gene to protein map based on
following criteria
1. if the mapping is to only one protein -- that is automatically taken
2. if the mapping is to more than one protein 
    a. it first consiers the canonical forms
    b. if there are multiple canonical forms matches it looks for psm
    c. if there are same psm matches, it looks for average intensity and chooses the protein that has maximum average intensity
    d. similarly it looks at isoforms and chooses the isoform with lowest number
    e. if multiple isoforms arise with same low number arises it does step 2b, 2c for isoforms too'''

def value_select(select_val_sp, psm_list,intensity_list):

    canonical = []
    canonical_index = []
    isoforms = []
    iso_index = []

    for i,new_val in enumerate(select_val_sp):
        if ("-" not in new_val) and ("sp|" in new_val):
            canonical.append(new_val)
            canonical_index.append(i)
        else:
            isoforms.append(new_val)
            iso_index.append(i)
#     print (isoforms)
    if len(canonical) == 1:
        final_value = canonical[0]
    if len(canonical) > 1:
        new_psm_list_canonical = [psm_list[index] for index in canonical_index]
        max_psm = max(new_psm_list_canonical)
        cnt_psm = new_psm_list_canonical.count(max_psm)

        if cnt_psm == 1:
            max_psm_index = new_psm_list_canonical.index(max_psm)
            final_value = select_val_sp[max_psm_index]
        else:
            new_intensity_list_canonical = [intensity_list[index] for index in canonical_index]
            max_intensity_index = new_intensity_list_canonical.index(np.max(new_intensity_list_canonical))
            final_value = select_val_sp[max_intensity_index]

    if (len(canonical) == 0) & (len(isoforms)==1):
        final_value = isoforms[0]

    if (len(canonical) == 0) & (len(isoforms)>1):

        iso_number = []

        for isoValue in isoforms:
            if "tr|" not in isoValue:
                iso_num=isoValue.split("-")[1][0]
                iso_number.append(iso_num)

                lowest_iso = iso_number.index(min(iso_number))
                cnt_low_iso = iso_number.count(min(iso_number))


                if cnt_low_iso == 1:
                    final_value = isoforms[lowest_iso]
        else:
            new_psm_list_isoforms = [psm_list[index] for index in iso_index]
            max_psm_iso = np.max(new_psm_list_isoforms)
            cnt_psm_iso = new_psm_list_isoforms.count(max_psm_iso)
            if cnt_psm_iso == 1:
                max_psm_iso_index = new_psm_list_isoforms.index(max_psm_iso)
                final_value = select_val_sp[max_psm_iso_index]
            else:
                new_intensity_list_isoforms = [intensity_list[index] for index in iso_index]
                max_intensity_iso_index = new_intensity_list_isoforms.index(np.max(new_intensity_list_isoforms))
                final_value = select_val_sp[max_intensity_iso_index] 
    #print (final_value)    
    return final_value



def select_iso(access_list, num_of_batches, intensity_list_avg, psm_list_avg):
    final_value = ""
    max_batches = max(num_of_batches)
    cnt_batches = num_of_batches.count(max_batches)
    if cnt_batches == 1:
        max_batch_index = num_of_batches.index(max_batches)
        final_value = access_list[max_batch_index]
#         print (final_value)
    else:
        batch_index_list =get_repeat_index(num_of_batches, max_batches)
#         print (batch_index_list)
        access_list_new, intensity_list_new, psm_list_new = select_correct_iso_lists(access_list, batch_index_list, intensity_list_avg, psm_list_avg)
        
#         print (access_list_new, intensity_list_new, psm_list_new)
        final_value = value_select(access_list_new, psm_list_new,intensity_list_new)
    if final_value == "":
        print (access_list, num_of_batches, intensity_list_avg, psm_list_avg)
    return final_value

        


def select_final_valtest_final(row):
    final_value = ""
    
    num_of_batches = list(row["nbatch"])
    intensity_list_avg = list(row["average_Intensity"])
    psm_list_avg = list(row["average_PSMs"])
    access_list =list(row["Protein Accession #"])
    
    if len(access_list) == 1:
        final_value = access_list[0]
    else:
        final_value = select_iso(access_list, num_of_batches, intensity_list_avg, psm_list_avg)
    return final_value


