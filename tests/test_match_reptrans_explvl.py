"""Tests for match representative transcript with expression level"""
import pytest
import pandas as pd
import numpy as np
from pandas.testing import assert_frame_equal
import tests.test_Functions as tFun
from transcript_sampler.match_reptrans_explvl import MatchReptransExplvl as match


class TestMatchReptrans:
    """Tests for match_reptrans_explvl.py"""
    # def test_gtf_to_df(self):
    # TO DO

    def test_dict_repr_trans_to_df(self):
        """
        This function test if a dict of {gene: representativeTranscript}
        is converted in a dataframe in the right format
        """
        dict_repr_test = {
            "ENSMUSG00000079415": "ENSMUST00000112933",
            "ENSMUSG00000024691": "ENSMUST00000025595",
            "ENSMUSG00000063683": "ENSMUST00000119960"}
        dict_mixed = {"a": 2, "b": 3}
        str_random = "jflkajflkaelfha"
        dict_int = {12: 34, 13: 66}
        data_frame = match.dict_repr_trans_to_df(dict_repr_test)
        datatype = {'Gene': np.dtype('O'), 'reprTrans': np.dtype('O')}

        with pytest.raises(TypeError, match=r"Only dictionaries are allowed"):
            match.dict_repr_trans_to_df(str_random)
        with pytest.raises(TypeError, match=r"Keys should be strings"):
            match.dict_repr_trans_to_df(dict_int)
        with pytest.raises(TypeError, match=r"Values should be strings"):
            match.dict_repr_trans_to_df(dict_mixed)

        assert tFun.column_number(data_frame) == 2, \
            "number of columns not equal to 2"
        assert tFun.column_dType(data_frame) == datatype, \
            "at least one column has the wrong datatype"
        assert tFun.duplicated_rows(data_frame).empty, \
            "at least one row is duplicated"
        assert tFun.NA_value(data_frame) == 0, \
            "at least one row contain NA values"

    def test_tsv_or_csv_to_df(self):
        """
        This function test if the function tsv_or_csv_to_df() can take
        csv and tsv file as input and return a pandas dataframe in the
        right format
        """
        path_tsv = tFun.find_path(r"test_gene_exprL")
        df_tsv = match.tsv_or_csv_to_df(path_tsv)
        path_csv = tFun.find_path(r"test_gene_exprL_csv.csv")
        df_csv = match.tsv_or_csv_to_df(path_csv)
        datatype = {'Transcript': np.dtype('O'),
                    'Expression_level': np.dtype('float64')}

        assert tFun.column_number(df_tsv) == 2, \
            "number of columns is not equal to 2"
        assert tFun.column_dType(df_tsv) == datatype, \
            "at least one column has the wrong datatype"
        assert tFun.duplicated_rows(df_tsv).empty, \
            "at least one row are duplicated "
        assert tFun.NA_value(df_tsv) == 0, \
            "at least one row contain NA values"
        assert_frame_equal(df_tsv, df_csv), \
            "csv and tsv import doesn't match"

    def test_expr_level_by_gene(self):
        """
        This function test if the function expr_level_by_gene can find
        the gene of each transcipt given by the expression level csv/tsv
        file and sum their expression level
        """
        path_tsv = tFun.find_path(r"test_gene_exprL")
        df_tsv_exprL = match.tsv_or_csv_to_df(path_tsv)
        df_gene_transcript = pd.DataFrame(
            {'Gene': ['ENSMUSG00000024691', 'ENSMUSG00000024691',
                      'ENSMUSG00000024691', 'ENSMUSG00000024691',
                      'ENSMUSG00000079415', 'ENSMUSG00000063683',
                      'ENSMUSG00000063683', 'ENSMUSG00000063683',
                      'ENSMUSG00000063683', 'ENSMUSG00000063683'],
             'Transcript': ['ENSMUST00000139270', 'ENSMUST00000151307',
                            'ENSMUST00000144662', 'ENSMUST00000025595',
                            'ENSMUST00000112933', 'ENSMUST000000449762',
                            'ENSMUST00000155846', 'ENSMUST00000157069',
                            'ENSMUST00000119960', 'ENSMUST00000123173']}
        )

        df_exprLevel = match.expr_level_by_gene(
            df_tsv_exprL, df_gene_transcript
            )
        datatype = {'Gene': np.dtype('O'),
                    'Expression_level': np.dtype('float64')}

        assert tFun.column_number(df_exprLevel) == 2, \
            "number of columns is not equal to 2"
        assert tFun.column_dType(df_exprLevel) == datatype, \
            "at least one column has the wrong datatype"
        assert tFun.duplicated_rows(df_exprLevel).empty, \
            "at least one row are duplicated "
        assert tFun.NA_value(df_exprLevel) == 0, \
            "at least one row contain NA values "
        assert tFun.duplicated_index(df_exprLevel).empty, \
            "at least one index element is duplicated"

    def test_match_by_gene(self):
        """
        This function test if the function "match_by_gene()" can
        create a pandas dataframe matching representative transcript
        and their expression level based on their gene in the
        correct pandas dataframe format.
        """

        dict_repr_test = {
            'ENSMUSG00000079415': 'ENSMUST00000112933',
            'ENSMUSG00000024691': 'ENSMUST00000025595',
            'ENSMUSG00000063683': 'ENSMUST00000119960'}
        df_dict_reprTrans = match.dict_repr_trans_to_df(dict_repr_test)

        path_tsv = tFun.find_path(r"test_gene_exprL")
        df_tsv_exprL = match.tsv_or_csv_to_df(path_tsv)
        df_gene_transcript = pd.DataFrame(
            {'Gene': ['ENSMUSG00000024691', 'ENSMUSG00000024691',
                      'ENSMUSG00000024691', 'ENSMUSG00000024691',
                      'ENSMUSG00000079415', 'ENSMUSG00000063683',
                      'ENSMUSG00000063683', 'ENSMUSG00000063683',
                      'ENSMUSG00000063683', 'ENSMUSG00000063683'],
             'Transcript': ['ENSMUST00000139270', 'ENSMUST00000151307',
                            'ENSMUST00000144662', 'ENSMUST00000025595',
                            'ENSMUST00000112933', 'ENSMUST000000449762',
                            'ENSMUST00000155846', 'ENSMUST00000157069',
                            'ENSMUST00000119960', 'ENSMUST00000123173']}
        )
        df_exprLevel = match.expr_level_by_gene(
            df_tsv_exprL, df_gene_transcript)

        df_match = match.match_by_gene(df_dict_reprTrans, df_exprLevel)
        datatype = {
            'reprTrans': np.dtype('O'),
            'Expression_level': np.dtype('float64')}

        assert tFun.column_number(df_match) == 2, \
            "number of columns is not equal to 2"
        assert tFun.column_dType(df_match) == datatype, \
            "at least one column has the wrong datatype"
        assert tFun.duplicated_rows(df_match).empty, \
            "at least one row are duplicated "
        assert tFun.NA_value(df_match) == 0, \
            "at least one row contain NA values "
        assert tFun.duplicated_index(df_match).empty, \
            "at least one index element is duplicated"

    def test_match_repr_transcript_expression_level(self):
        """
        This function test that the right output is generated by the function
        match_repr_transcript_expression_level()
        """
        input_path = tFun.find_path("test_gene_exprL")
        intermediate_path = tFun.find_path_intermediateFile()
        dict_repr_test = {
            'ENSMUSG00000079415': 'ENSMUST00000112933',
            "ENSMUSG00000024691": "ENSMUST00000025595",
            "ENSMUSG00000063683": "ENSMUST00000119960"}

        match.match_repr_transcript_expression_level(
            exprTrans=input_path,
            dict_reprTrans=dict_repr_test,
            gtf_file=intermediate_path)

        ref_path = tFun.find_path("test_ref_output.tsv")
        output_path = tFun.find_output()

        with open(ref_path, 'r', encoding="utf-8") as t1,\
            open(output_path, 'r', encoding="utf-8") as t2,\
            open(input_path, 'r', encoding="utf-8") as t3:
            fileRef = t1.readlines()
            fileOutput = t2.readlines()
            fileInput = t3.readlines()

        assert (
            sorted(fileRef) == sorted(fileOutput)
            ), "the output does't match the expected tsv file"
        assert (
            sorted(fileRef) != sorted(fileInput)
            ), "the output does't match the expected tsv file"

    # def test_txt_to_dict(self):
    #     """This function tests if txt is convertod to dict"""
    #     path = tFun.find_path("test_dict_repr_trans.txt")
    #     dico = match.txt_to_dict(path)
    #     dict_test = {'ENSMUSG00000079415': 'ENSMUST00000112933',
    #                 "ENSMUSG00000024691": "ENSMUST00000025595",
    #                 "ENSMUSG00000063683": "ENSMUST00000119960"}
    #     assert dico == dict_test

    # def test_transcripts_by_gene_inDf():
    #     """
    #     This function test if a dataframe generated from
    #     the intermediate file is converted in another
    #     dataframe without the support level column.
    #     """
    #     path = tFun.find_path_intermediateFile()
    #     df = repr.import_gtfSelection_to_df(path)
    #     df_gene = match.transcripts_by_gene_inDf(df)
    #     datatype = {'Gene': np.dtype('O'), 'Transcript': np.dtype('O')}
    #     assert tFun.column_number(df_gene) == (
    #         2, "number of columns is not equal to 2")
    #     assert tFun.column_dType(df_gene) == (
    #         datatype, "at least one column has the wrong datatype")
    #     assert tFun.duplicated_rows(df_gene).empty, \
    #         "at least one row are duplicated"
    #     assert tFun.NA_value(df_gene) == 0, "at least one row contain NA values"

    # def test_output_tsv():
    #     """
    #     This function test if a tsv file is generated from a pandas
    #     dataframe in the right format.
    #     """

    #     dict_repr_test = {
    #         'ENSMUSG00000079415': 'ENSMUST00000112933',
    #         "ENSMUSG00000024691": "ENSMUST00000025595",
    #         "ENSMUSG00000063683": "ENSMUST00000119960"}
    #     df_dict_reprTrans = match.dict_repr_trans_to_df(dict_repr_test)

    #     path_tsv = tFun.find_path(r"test_gene_exprL")
    #     df_tsv_exprL = match.tsv_or_csv_to_df(path_tsv)
    #     path_intermediate = tFun.find_path_intermediateFile()
    #     df_intermediate = repr.import_gtfSelection_to_df(path_intermediate)
    #     df_gene_transcript = match.transcripts_by_gene_inDf(df_intermediate)

    #     df_exprLevel = match.expr_level_by_gene(df_tsv_exprL, df_gene_transcript)

    #     df_match = match.match_by_gene(df_dict_reprTrans, df_exprLevel)

    #     match.output_tsv(df_match)

    #     ref_path = tFun.find_path("test_ref_output.tsv")
    #     output_path = tFun.find_output()

    #     with open(ref_path, 'r') as t1, open(output_path, 'r') as t2:
    #         fileRef = t1.readlines()
    #         fileOutput = t2.readlines()

    #     assert (
    #         sorted(fileRef) == sorted(fileOutput)
    #         ), "the output does't match the expected tsv file"

# test_dict_repr_trans_to_df()
# test_txt_to_dict()
# test_transcripts_by_gene_inDf()
# test_tsv_or_csv_to_df()
# test_expr_level_by_gene()
# test_match_by_gene()
# test_output_tsv()
# test_match_repr_transcript_expression_level()

# print("test_match is done ! No error was found")
