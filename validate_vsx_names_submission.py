import re
import pandas as pd


def validate_vsx_names_submission(csv_path="data/others/vsx_names_submissions_draft.csv"):
    """Validate the format of the VSX Names submission spreadsheet."""

    def validate_other_names(names_str):
        # the names need to be comma separated, without spaces  between them
        names = names_str.split(",")
        for name in names:
            if re.search(r"^\s+|\s$", name):
                return False  # leading / trailing spaces. invalid
        return True

    def is_non_empty_str(a_str):
        return not pd.isna(a_str) and a_str != ""

    def print_df_with_defaults(df):
        with pd.option_context("display.max_rows", None, "display.expand_frame_repr", False, "display.max_columns", None):
            print(df[["VSX Name", "Other names", "oid"]])

    #
    # the main logic
    df = pd.read_csv(csv_path)

    df_invalid_other_names = df[~pd.Series([validate_other_names(n) for n in df["Other names"]])]
    if len(df_invalid_other_names) > 0:
        print("Other names in some entries are invalid:")
        print_df_with_defaults(df_invalid_other_names)
        return None

    df_missing_vsx_name = df[~pd.Series([is_non_empty_str(n) for n in df["VSX Name"]])]
    if len(df_missing_vsx_name) > 0:
        print("VSX Name in some entries are missing:")
        print_df_with_defaults(df_missing_vsx_name)
        return None

    #  case all okay
    print(f"All {len(df)} entries valid")
    return df


df = validate_vsx_names_submission()
df
