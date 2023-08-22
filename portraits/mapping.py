import pandas as pd
import logging

def get_gs_for_probes_from_3col(platform_file, probe_list):
    
    import logging

    """
    Getting probe-gene symbol dictionary

    :param platform_name: str, platform name
    :param probe_list: list, list with probe names

    :return: dict, dictionary with probe-gene symbol key-values
    """
    try:
        platform_data = pd.read_csv(platform_file, sep='\t', header=None, index_col=0, na_values=["NONE"])
    except Exception as e:
        logging.warning(f"Failed to read mapping 3col-file: {str(e)}")
        return None

    dict_raw_name_id = dict()
    not_found_probes_amount = len(set(probe_list).difference(platform_data.index))
    if not_found_probes_amount:
        logging.warn(f'{not_found_probes_amount} probes not found or format is not correct.')
        return dict()

    result = platform_data[1].loc[probe_list].dropna().astype(str).apply(
        lambda x: x.strip().replace(" ", "").split("///")).to_dict()

    return result


def get_expressions_list(probes_list, probes_value_table, method='max'):
    """
    Returns list of expressions for matching gene_symbol

    :param probes_list: list, list of probes (for matching gene_symbol)
    :param probes_value_table: pd.DataFrame, matching table for probe_id and expression values (for each sample)
    :param method: str, getting expressions method (max / med )

    :return: list, list of expressions for matching gene_symbol
    """

    def average_expression(expressions_list):
        return sum(expressions_list) / len(expressions_list)

    probes_avg_expr_dict = {}
    # count average for all gsms
    for probe in probes_list:
        probes_avg_expr_dict[probe] = average_expression(probes_value_table.loc[probe, :])
    probe_res = ''
    if method == 'max':
        # choose probe-id with max average value
        probe_res = max(probes_avg_expr_dict, key=probes_avg_expr_dict.get)

    elif method == 'med':
        # choose probe with median value
        # if 2 probes choose probe with max average value

        # sort dict by values and return list of keys
        sorted_probes_list = sorted(probes_avg_expr_dict, key=probes_avg_expr_dict.get)
        if len(probes_list) % 2 == 0:
            probe_res = max(sorted_probes_list[len(probes_list) / 2 - 1],
                            sorted_probes_list[len(probes_list) / 2],
                            key=probes_avg_expr_dict.get)
        else:
            probe_res = sorted_probes_list[len(probes_list) / 2]
    return probes_value_table.loc[probe_res, :]

    
def get_expressions_for_gs(probes_gs_dict, probes_value_table, gs_sel_alg='max'):
    """
    Getting genes/samples expression table

    :param probes_gs_dict: dict, dictionary with probe-gene symbol key-values
    :param probes_value_table: pd.DataFrame, probes/samples expression transformed dataframe
    :param gs_sel_alg: str, getting expressions method (max / med )

    :return: gs_expr_table: pd.DataFrame, genes/samples expression table
    """
    import logging
    
    def get_reverse_dictionary(probes_gs_dict):
        from collections import defaultdict
        gs_probes_dict = defaultdict(list)
        for probe, gs_list in probes_gs_dict.items():
            for gs in gs_list:
                gs_probes_dict[gs].append(probe)
        return gs_probes_dict
    
    gs_expr_table = pd.DataFrame()

    logging.info("Making list of probes for each of gene-symbols ...")
    gs_probes_dict = get_reverse_dictionary(probes_gs_dict)

    logging.info("Making expression list for each of gene-symbols ...")
    for gs, probe_list in gs_probes_dict.items():
        # print(probe_list[:30])
        if (len(probe_list) == 1):
            gs_expr_table[gs] = probes_value_table.loc[probe_list[0], :]
        else:
            gs_expr_table[gs] = get_expressions_list(
                probes_list=probe_list,
                probes_value_table=probes_value_table,
                method=gs_sel_alg
            )

    logging.info(f'Final expression samples/gene symbols table shape: {gs_expr_table.shape}')

    return gs_expr_table