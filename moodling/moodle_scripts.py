def combine_fb_and_report(cand_num, orig_name, fb_name):
    mergedObject = PdfFileMerger()
    mergedObject.append(PdfFileReader(fb_dir + '/' + cand_num+ '.pdf', 'rb'))
    mergedObject.append(PdfFileReader(report_dir + '/' + orig_name, 'rb'))
    mergedObject.write(dst_dir + '/' + fb_name)
    return(fb_name)
    
    
def create_feedback_file(anon_id, orig_fn, cand_num):
    cn_name = cand_num + '.pdf'
    cn_path = os.path.join(sh_dir, cn_name)

    orig_fn = orig_fn[:-5]
    new_name = 'Participant_' + anon_id + '_assignsubmission_file_' + orig_fn + '_feedback.pdf'
    new_path = os.path.join(dst_dir, new_name)
    
    #print(cn_path, new_path)
    #print(fb_file_path)
    shutil.copy(cn_path, new_path)
    return(new_name)
    
def get_cand_num(orig_fn):
    match = re.search('\d{9}', orig_fn)
    if match == None:
        match = re.search('\d{5}', orig_fn)
    cn = match.group(0)
    if len(cn)==5:
        cn = 'CN' + cn
    else:
        cn = 'SN' + cn
    return(cn)
    
def get_anon_id(fn):
    anon_id = fn.split('_')[1]
    return(anon_id)