import numpy as np
from astropy.table.pprint import TableFormatter


# * Special table formatting for bounded (val, min, max) values
def fmt_func(fmt):
    def _fmt(x):
        ret = f"{x['val']:{fmt}}"
        if np.isnan(x['min']) and np.isnan(x['max']):
            return ret + " (fixed)"
        else:
            mn = ("-∞" if np.isnan(x['min']) or x['min'] == -np.inf
                  else f"{x['min']:{fmt}}")
            mx = ("∞" if np.isnan(x['max']) or x['max'] == np.inf
                  else f"{x['max']:{fmt}}")
            return f"{ret} ({mn}, {mx})"
    return _fmt


class BoundedParTableFormatter(TableFormatter):
    """Format bounded parameters.
    Bounded parameters are 3-field structured arrays, with fields
    'val', 'min', and 'max'.  To be set as Table's `TableFormatter'.
    """
    def _pformat_table(self, table, *args, **kwargs):
        bpcols = []
        try:
            for col in table.columns.values():
                if len(col.dtype) == 3:  # bounded!
                    bpcols.append((col, col.info.format))
                    col.info.format = fmt_func(col.info.format or "g")
            return super()._pformat_table(table, *args, **kwargs)
        finally:
            for col, fmt in bpcols:
                col.info.format = fmt
