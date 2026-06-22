#After running PAHFIT model.fit you can use 'model' as the input to pahfit_to_ecsv to save all fit params with associated errors from covariance 
#matrices. The rest of the functions in here are just helper functions to make pahfit_to_ecsv work. -Hannah

import numpy as np
from astropy.table import Table


def build_flat_name_map(astropy_model):
    return {v: k for k, v in astropy_model._param_map.items()}
def get_covariance(model, comp1, par1, comp2, par2):
    astropy_model = model.fitter.model

    flat_name_map = build_flat_name_map(astropy_model)

    # component index
    idx1 = astropy_model.submodel_names.index(comp1)
    idx2 = astropy_model.submodel_names.index(comp2)

    # map to flat names
    flat1 = flat_name_map[(idx1, par1)]
    flat2 = flat_name_map[(idx2, par2)]

    return astropy_model.cov_matrix[flat1, flat2]
def get_power_and_error(model, feature_name):
    astropy_model = model.fitter.model

    comp = astropy_model[feature_name]
    power_val = comp.power.value

    try:
        var = get_covariance(model, feature_name, "power",
                                   feature_name, "power")
        power_err = np.sqrt(var)
    except KeyError:
        power_err = None

    return {"power": power_val, "power_err": power_err}

def pahfit_to_ecsv(model,
                  output_file=None,
                  include_components=None,
                  exclude_components=None,
                  include_params=None,
                  add_metadata=True,
                  verbose=False):
    """
    Extract PAHFIT model parameters + uncertainties from covariance matrix
    and optionally save to an ECSV file.

    Parameters
    ----------
    model : PAHFIT model object
        The fitted PAHFIT model.

    output_file : str, optional
        Path to save ECSV file. If None, does not save.

    include_components : list of str, optional
        Only include these component names.

    exclude_components : list of str, optional
        Skip these component names.

    include_params : list of str, optional
        Only include these parameter names.

    add_metadata : bool
        Whether to include fit metadata (chi2, etc.)

    verbose : bool
        Print progress.

    Returns
    -------
    table : astropy.table.Table
    """

    astropy_model = model.fitter.model

    # ---- helper: map flattened covariance indices ----
    def build_flat_name_map(astropy_model):
        return {v: k for k, v in astropy_model._param_map.items()}

    flat_name_map = build_flat_name_map(astropy_model)

    def get_covariance(comp1, par1, comp2, par2):
        try:
            idx1 = astropy_model.submodel_names.index(comp1)
            idx2 = astropy_model.submodel_names.index(comp2)

            flat1 = flat_name_map[(idx1, par1)]
            flat2 = flat_name_map[(idx2, par2)]

            return astropy_model.cov_matrix[flat1, flat2]
        except Exception:
            return np.nan

    def get_param_and_error(comp_name, param_name):
        comp = astropy_model[comp_name]

        try:
            value = getattr(comp, param_name).value
        except AttributeError:
            return np.nan, np.nan

        var = get_covariance(comp_name, param_name,
                             comp_name, param_name)

        if var is None or np.isnan(var) or var < 0:
            err = np.nan
        else:
            err = np.sqrt(var)

        return value, err

    # ---- build table ----
    rows = []

    for comp_name in astropy_model.submodel_names:

        # Apply filters
        if include_components and comp_name not in include_components:
            continue
        if exclude_components and comp_name in exclude_components:
            continue

        comp = astropy_model[comp_name]

        if not hasattr(comp, "param_names"):
            continue

        if verbose:
            print(f"Processing {comp_name}")

        row = {"feature": comp_name}

        for param in comp.param_names:

            if include_params and param not in include_params:
                continue

            val, err = get_param_and_error(comp_name, param)

            row[param] = val
            row[f"{param}_err"] = err

        rows.append(row)

    table = Table(rows)

    # ---- add metadata ----
    if add_metadata:
        try:
            table.meta["chi2"] = model.fit_info.get("chi2", np.nan)
            table.meta["redchi"] = model.fit_info.get("redchi", np.nan)
        except Exception:
            pass

    # ---- save if requested ----
    if output_file is not None:
        table.write(output_file, format="ascii.ecsv", overwrite=True)
        if verbose:
            print(f"Saved to {output_file}")

    return table





def covariance_to_table(model):
    astropy_model = model.fitter.model
    cov = astropy_model.cov_matrix

    names = list(getattr(cov, "param_names", astropy_model.param_names))
    cov_array = np.array(cov.cov_matrix if hasattr(cov, "cov_matrix") else cov, dtype=float)

    table = Table(cov_array, names=names)
    table.add_column(names, name="parameter", index=0)

    return table

def features_with_uncertainties(model, precision=6, include_covariance_meta=False):
    param_unc_table = pahfit_to_ecsv(model)
    final_features = model.features.copy()

    for col in ["temperature_unc", "tau_unc", "power_unc", "wavelength_unc", "fwhm_unc"]:
        if col not in final_features.colnames:
            final_features[col] = np.full(len(final_features), np.nan)
        else:
            final_features[col][:] = np.nan

    for col in ["temperature_pm", "tau_pm", "power_pm", "wavelength_pm", "fwhm_pm"]:
        if col not in final_features.colnames:
            final_features[col] = np.full(len(final_features), "", dtype="U64")
        else:
            final_features[col][:] = ""

    def param_value(row, name):
        try:
            return float(row[name])
        except Exception:
            return np.nan

    def pm_string(val, err):
        if not np.isfinite(val):
            return ""
        if np.isfinite(err):
            return f"{val:.{precision}g} ± {err:.{precision}g}"
        return f"{val:.{precision}g}"

    def assign(idx, pname, val, err):
        final_features[f"{pname}_unc"][idx] = err
        final_features[f"{pname}_pm"][idx] = pm_string(val, err)

    for row in param_unc_table:
        feature_name = str(row["feature"])
        idx = np.where(final_features["name"] == feature_name)[0]
        if len(idx) == 0:
            continue

        kind = str(final_features["kind"][idx[0]])

        if kind in ("starlight", "dust_continuum"):
            assign(idx, "temperature", param_value(row, "temperature"), param_value(row, "temperature_err"))
            assign(idx, "tau", param_value(row, "amplitude"), param_value(row, "amplitude_err"))

        elif kind == "line":
            assign(idx, "power", param_value(row, "power"), param_value(row, "power_err"))
            assign(idx, "wavelength", param_value(row, "mean"), param_value(row, "mean_err"))
            assign(idx, "fwhm", 2.355 * param_value(row, "stddev"), 2.355 * param_value(row, "stddev_err"))

        elif kind == "dust_feature":
            assign(idx, "power", param_value(row, "power"), param_value(row, "power_err"))
            assign(idx, "wavelength", param_value(row, "x_0"), param_value(row, "x_0_err"))
            assign(idx, "fwhm", param_value(row, "fwhm"), param_value(row, "fwhm_err"))

        elif kind == "attenuation":
            assign(idx, "tau", param_value(row, "tau_sil"), param_value(row, "tau_sil_err"))

        elif kind == "absorption":
            assign(idx, "tau", param_value(row, "tau"), param_value(row, "tau_err"))
            assign(idx, "wavelength", param_value(row, "x_0"), param_value(row, "x_0_err"))
            assign(idx, "fwhm", param_value(row, "fwhm"), param_value(row, "fwhm_err"))

    if include_covariance_meta:
        cov = model.fitter.model.cov_matrix
        final_features.meta["covariance_parameter_order"] = list(getattr(cov, "param_names", model.fitter.model.param_names))
        final_features.meta["covariance_matrix"] = np.array(cov.cov_matrix if hasattr(cov, "cov_matrix") else cov, dtype=float).tolist()

    return final_features