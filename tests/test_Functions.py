import pandas as pd
import numpy as np
import os


def find_path(filename: str) -> str:
    """Find the path to a file

        Args:
            name of a file

        Returns:
            str path of a file

        Raises:
            None
    """
    absolute_path = os.path.dirname(__file__)
    test_file = "inputs/" + str(filename)
    full_path = os.path.join(absolute_path, test_file)
    return full_path


def find_output():
    """Find the path of the output file

        Args:
            name of a file

        Returns:
            str path of a file

        Raises:
            None
    """
    absolute_path = os.path.dirname(__file__)
    test_file = "ReprTrans_ExpressionLevel.tsv"
    full_path = os.path.join(absolute_path, test_file)
    return full_path


def find_path_intermediateFile() -> str:
    """Find the path to gencode.vM31.annotation_intermediat_file.txt

        Args:
            none

        Returns:
            str path of gencode.vM31.annotation_intermediat_file.txt

        Raises:
            None
    """
    absolute_path = os.path.dirname(__file__)
    test_file = r"inputs/test_gencode.vM31.annotation_intermediat_file.txt"
    full_path = os.path.join(absolute_path, test_file)
    return full_path


def column_number(df: pd.DataFrame) -> int:

    """Return the number of column of a df

        Args:
            dataframe

        Returns:
            int

        Raises:
            None
    """
    length = len(df.columns)
    return length


def column_dType(df: pd.DataFrame) -> dict[str, np.dtype]:
    """Return the type of each column of a df in a dict

        Args:
            Pandas dataframe

        Returns:
            dict{column:np.dtype()}

        Raises:
            None
    """
    dtype = df.dtypes.to_dict()
    return dtype


def duplicated_rows(df: pd.DataFrame) -> pd.DataFrame:
    """Return the sum of duplicated rows in a df

        Args:
            Pandas dataframe

        Returns:
            int

        Raises:
            None
    """
    df_dupl = df[df.duplicated()]
    return df_dupl


def duplicated_index(df: pd.DataFrame) -> pd.DataFrame:
    """Return the sum of duplicated index in a df

        Args:
            Pandas dataframe

        Returns:
            int

        Raises:
            None
    """
    df_dupl = df[df.index.duplicated()]
    return df_dupl


def NA_value(df: pd.DataFrame) -> int:
    """Return the sum of NA values in a df

        Args:
            Pandas dataframe

        Returns:
            int

        Raises:
            None
    """
    nNA = df.isna().sum().sum()
    return nNA
