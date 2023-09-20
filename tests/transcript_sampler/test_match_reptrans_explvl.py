"""Tests for match representative transcript with expression level."""
import pytest  # type: ignore
import pandas as pd  # type: ignore
import numpy as np  # type: ignore
from pandas.testing import assert_frame_equal  # type: ignore

from ...scRNAsim_toolz.transcript_sampler.match_explvl import (  # type: ignore
    MatchReptransExplvl as match
)
import test_functions as tFun  # type: ignore


class TestMatchReptrans:
    """Tests for match_reptrans_explvl.py."""

    def test_dict_repr_trans_to_df(self):
        """Test dict_repr_trans_to_df() function.

        This function test if a dict of {gene: representativeTranscript}.
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
        assert tFun.column_d_type(data_frame) == datatype, \
            "at least one column has the wrong datatype"
        assert tFun.duplicated_rows(data_frame).empty, \
            "at least one row is duplicated"
        assert tFun.na_value(data_frame) == 0, \
            "at least one row contains NA values"

    def test_tsv_or_csv_to_df(self):
        """Test tsv_or_csv_to_df() function.

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
        assert tFun.column_d_type(df_tsv) == datatype, \
            "at least one column has the wrong datatype"
        assert tFun.duplicated_rows(df_tsv).empty, \
            "at least one row is duplicated"
        assert tFun.na_value(df_tsv) == 0, \
            "at least one row contains NA values"
        assert assert_frame_equal(df_tsv, df_csv) is None, \
            "csv and tsv import doesn't match"

    def test_expr_level_by_gene(self):
        """Test expr_level_by_gene() function.

        This function test if the function expr_level_by_gene can find
        the gene of each transcript given by the expression level csv/tsv
        file and sum their expression level
        """
        path_tsv = tFun.find_path(r"test_gene_exprL")
        df_tsv_exp_lvl = match.tsv_or_csv_to_df(path_tsv)
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

        df_exp_lvl = match.expr_level_by_gene(
            df_tsv_exp_lvl, df_gene_transcript
            )
        datatype = {'Gene': np.dtype('O'),
                    'Expression_level': np.dtype('float64')}

        assert tFun.column_number(df_exp_lvl) == 2, \
            "number of columns is not equal to 2"
        assert tFun.column_d_type(df_exp_lvl) == datatype, \
            "at least one column has the wrong datatype"
        assert tFun.duplicated_rows(df_exp_lvl).empty, \
            "at least one row is duplicated"
        assert tFun.na_value(df_exp_lvl) == 0, \
            "at least one row contains NA values"
        assert tFun.duplicated_index(df_exp_lvl).empty, \
            "at least one index element is duplicated"

    def test_match_by_gene(self):
        """Test match_by_gene() function.

        This function test if the function "match_by_gene()" can
        create a pandas dataframe matching representative transcript
        and their expression level based on their gene in the
        correct pandas dataframe format.
        """
        dict_repr_test = {
            'ENSMUSG00000079415': 'ENSMUST00000112933',
            'ENSMUSG00000024691': 'ENSMUST00000025595',
            'ENSMUSG00000063683': 'ENSMUST00000119960'}
        df_dict_repr_trans = match.dict_repr_trans_to_df(dict_repr_test)

        path_tsv = tFun.find_path(r"test_gene_exprL")
        df_tsv_exp_lvl = match.tsv_or_csv_to_df(path_tsv)
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
        df_exp_lvl = match.expr_level_by_gene(
            df_tsv_exp_lvl, df_gene_transcript)

        df_match = match.match_by_gene(df_dict_repr_trans, df_exp_lvl)
        datatype = {
            'reprTrans': np.dtype('O'),
            'Expression_level': np.dtype('float64')}

        assert tFun.column_number(df_match) == 2, \
            "number of columns is not equal to 2"
        assert tFun.column_d_type(df_match) == datatype, \
            "at least one column has the wrong datatype"
        assert tFun.duplicated_rows(df_match).empty, \
            "at least one row is duplicated"
        assert tFun.na_value(df_match) == 0, \
            "at least one row contains NA values"
        assert tFun.duplicated_index(df_match).empty, \
            "at least one index element is duplicated"

    def test_match_repr_transcript_expression_level(self):
        """Test match_repr_transcript_expression_level().

        This function test that the right output is generated by the
        function match_repr_transcript_expression_level().
        """
        input_path = tFun.find_path("test_gene_exprL")
        gtf_file = tFun.find_path("test.gtf")
        dict_repr_test = {
            'ENSMUSG00000079415': 'ENSMUST00000112933',
            "ENSMUSG00000024691": "ENSMUST00000025595",
            "ENSMUSG00000063683": "ENSMUST00000119960"}

        # Create an instance of MatchReptransExplvl
        match_instance = match()

        df_result = match_instance.match_repr_transcript_expression_level(
            expr_trans=input_path,
            dict_repr_trans=dict_repr_test,
            gtf_file=gtf_file
        )

        ref_path = tFun.find_path("test_ref_output.tsv")
        output_path = tFun.find_output()

        with open(
            ref_path, 'r', encoding="utf-8"
        ) as test_file_1, open(
            output_path, 'r', encoding="utf-8"
        ) as test_file_2:
            file_ref = test_file_1.readlines()
            file_output = test_file_2.readlines()

        assert (
            sorted(file_ref) == sorted(file_output)
        ), "the output doesn't match the expected tsv file"
        assert (
            sorted(file_ref) != sorted(
                df_result.to_csv(index=False).splitlines()
            )
        ), "the output doesn't match the expected tsv file"
