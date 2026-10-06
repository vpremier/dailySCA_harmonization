#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 11:08:04 2026

@author: vpremier
"""

import pandas as pd
import numpy as np







def estimate_snow_cover_status(
    df,
    epsilon=10.0,
    alpha=0.5  # interpolation strength
):
    """
    Estimate snow cover using status_ts as physical constraint.
    status_ts must be in same units as SCA (e.g. % of basin).
    """

    df = df.copy()

    S_est = []
    trends = []

    S_prev = None

    for i, row in df.iterrows():
        Smin = row["sca_min"]
        Smax = row["sca_max"]
        S_snow = row["status_ts"]  # assumed already normalized
        # if i==1:
            
        #     dd
        width = Smax - Smin
        
        # --- Case 1: trustworthy ---
        if width <= epsilon:
            S = (Smin + Smax) / 2.0
            trend = "trusted"

        else:
            if S_prev is None:
                S = (Smin + Smax) / 2.0
                trend = "init"

            else:
                # --- Determine trend from status ---
                if S_snow > S_prev:
                    trend = "up"
                elif S_snow < S_prev:
                    trend = "down"
                else:
                    trend = "flat"

                # --- Apply constraints ---
                if trend == "up":
                    upper = min(Smax, S_snow)
                    lower = max(Smin, S_prev)

                    if upper < lower:
                        # inconsistent → fallback
                        S = upper
                    else:
                        S = (1 - alpha) * lower + alpha * upper

                elif trend == "down":
                    upper = min(Smax, S_prev)
                    lower = Smin

                    if upper < lower:
                        S = lower
                    else:
                        S = (1 - alpha) * upper + alpha * lower

                else:
                    S = np.clip(S_prev, Smin, Smax)

        S_est.append(S)
        trends.append(trend)
        S_prev = S

    df["sca_est"] = S_est
    df["trend"] = trends

    return df



def estimate_snow_cover(
    df,
    epsilon=10.0,           # % threshold for trusting range
    Tsnow=10.0,            # snow threshold (°C)
    Tmelt=10.0,            # melt threshold (°C)
    precip_thresh=5,    # minimum effective precipitation
    alpha_up=0.1,         # weight toward upper bound when increasing
    alpha_down=0.1        # weight toward lower bound when decreasing
):
    """
    Estimate snow cover from min/max bounds using physical rules.
    Returns a new DataFrame with estimated SCA and trend.
    """

    df = df.copy()

    # Compute effective precipitation
    df["precip_eff"] = df["pr_mean"] * df["status_ts"]

    # Initialize outputs
    S_est = []
    trends = []

    S_prev = None

    for i, row in df.iterrows():
        Smin = row["sca_min"]
        Smax = row["sca_max"]
        T = row["ta_mean"] - 273.15
        P = row["pr_mean"]

        width = Smax - Smin

        # --- Case 1: trustworthy range ---
        if width <= epsilon:
            S = (Smin + Smax) / 2.0
            trend = "trusted"

        else:
            # --- Determine trend ---
            if (P > precip_thresh) and (T <= Tsnow):
                trend = "up"
            elif (P <= precip_thresh) and (T > Tmelt):
                trend = "down"
            else:
                trend = "flat"

            # --- First timestep fallback ---
            if S_prev is None:
                S = (Smin + Smax) / 2.0

            else:
                # --- Apply constraints ---
                if trend == "up":
                    lower = max(S_prev, Smin)
                    S = (1 - alpha_up) * lower + alpha_up * Smax

                elif trend == "down":
                    upper = min(S_prev, Smax)
                    S = alpha_down * Smin + (1 - alpha_down) * upper

                else:  # flat / uncertain
                    S = np.clip(S_prev, Smin, Smax)

        # Store
        S_est.append(S)
        trends.append(trend)
        S_prev = S

    df["sca_est"] = S_est
    df["trend"] = trends

    return df

import matplotlib.pyplot as plt

def plot_snow_cover(df):
    plt.figure(figsize=(10, 5))

    x = df.index

    # Shaded min-max range
    plt.fill_between(
        x,
        df["sca_min"],
        df["sca_max"],
        alpha=0.2,
        label="SCA range (min–max)"
    )

    # Estimated snow cover
    plt.plot(
        x,
        df["sca_est"],
        linewidth=2,
        label="Estimated SCA"
    )

    # Optional: midpoint for reference
    midpoint = (df["sca_min"] + df["sca_max"]) / 2
    plt.plot(
        x,
        midpoint,
        linestyle="--",
        linewidth=1,
        label="Midpoint (naive)"
    )

    plt.xlabel("Time")
    plt.ylabel("Snow Cover (%)")
    plt.title("Snow Cover Estimation from Uncertain Bounds")
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.show()
    


import matplotlib.pyplot as plt

def plot_all(df):
    fig, ax1 = plt.subplots(figsize=(11, 6))

    x = df.index

    # --- Snow cover (primary axis) ---
    ax1.fill_between(x, df["sca_min"], df["sca_max"], alpha=0.2, label="SCA range")
    ax1.plot(x, df["sca_est"], linewidth=2, label="Estimated SCA")
    ax1.set_ylabel("Snow Cover (%)")
    ax1.set_xlabel("Time")

    # --- Temperature (secondary axis) ---
    ax2 = ax1.twinx()
    ax2.plot(x, df["ta_mean"]-273.15, linestyle="--", label="Temperature")
    ax2.set_ylabel("Temperature (°C)")

    # --- Precipitation (third axis) ---
    ax3 = ax1.twinx()
    ax3.spines["right"].set_position(("outward", 60))  # shift axis outward
    ax3.bar(x, df["pr_mean"], alpha=0.7, label="Precipitation")
    ax3.set_ylabel("Precipitation")

    # --- Combine legends ---
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    lines_3, labels_3 = ax3.get_legend_handles_labels()

    ax1.legend(
        lines_1 + lines_2 + lines_3,
        labels_1 + labels_2 + labels_3,
        loc="upper right"
    )

    ax1.grid(True)
    plt.title("Snow Cover, Temperature, and Precipitation")

    plt.tight_layout()
    plt.show()



    
df1= pd.read_csv(r'/mnt/CEPH_PROJECTS/SNOWCOP/Vale/daily_test/Area06_hy2122_meteo.csv')
df2= pd.read_csv(r'/mnt/CEPH_PROJECTS/SNOWCOP/Vale/daily_test/Area06_hy2122_sca.csv')

df = pd.merge(df1, df2, on="time", how="inner")
df = df.fillna(0)


result = estimate_snow_cover_status(df)
plot_snow_cover(result)
plot_all(result)
