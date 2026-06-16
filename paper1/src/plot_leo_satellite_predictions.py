"""LEO 卫星通信质量和可视度预测绘图脚本。

根据文档中的模型和注释，生成图1~图10的绘制代码。
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import warnings

# 配置 Matplotlib 中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 压制字体相关的 UserWarning
warnings.filterwarnings('ignore', category=UserWarning, module='matplotlib')

# 常量（单位：km、s）
RE = 6371.0
h = 520.0
g = 9.81 / 1000.0  # km/s^2
Ts = np.sqrt(4 * np.pi**2 * (RE + h)**3 / (g * RE**2))
Tg = 86164.098903691

# 卫星轨道参数
inclination_deg = 96.5
omega = 0.0
lat_gs_deg = -33.92
lon_gs_deg = 18.86
hs = 0.111  # km

# 地球椭球参数
a_ell = 6378.137
b_ell = 6356.752
epsilon2 = (a_ell**2 - b_ell**2) / a_ell**2

# 时间设置
t_end = 3 * 86400.0
dt = 30.0
k_steps = int(t_end / dt) + 1

t = np.linspace(0, t_end, k_steps)

# 最大通信范围
d_max = RE * np.tan(np.arccos(RE / (RE + h)))


def rotation_matrix_x(angle_rad):
    c, s = np.cos(angle_rad), np.sin(angle_rad)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])


def rotation_matrix_y(angle_rad):
    c, s = np.cos(angle_rad), np.sin(angle_rad)
    return np.array([[c, 0.0, -s], [0.0, 1.0, 0.0], [s, 0.0, c]])


def rotation_matrix_z(angle_rad):
    c, s = np.cos(angle_rad), np.sin(angle_rad)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def ground_station_initial_position():
    lat = np.deg2rad(lat_gs_deg)
    lon = np.deg2rad(lon_gs_deg)
    Rtrans = a_ell / np.sqrt(1.0 - epsilon2 * np.sin(lat)**2)
    x = (Rtrans + hs) * np.cos(lat) * np.cos(lon)
    y = (Rtrans + hs) * np.cos(lat) * np.sin(lon)
    z = ((1 - epsilon2) * Rtrans + hs) * np.sin(lat)
    return np.array([x, y, z])


def satellite_positions():
    inclination = np.deg2rad(inclination_deg)
    s0 = np.array([0.0, 0.0, RE + h])
    L = rotation_matrix_y(inclination)
    H = rotation_matrix_z(omega)
    S_step = 2.0 * np.pi * dt / Ts
    Q = rotation_matrix_x(S_step)
    Qt = np.eye(3)
    positions = np.zeros((k_steps, 3))

    for k in range(k_steps):
        S = H @ L @ Qt @ s0
        positions[k] = S
        Qt = Q @ Qt

    return positions


def ground_station_positions():
    g0 = ground_station_initial_position()
    G_step = 2.0 * np.pi * dt / Tg
    R = rotation_matrix_z(G_step)
    Rt = np.eye(3)
    positions = np.zeros((k_steps, 3))

    for k in range(k_steps):
        positions[k] = Rt @ g0
        Rt = R @ Rt

    return positions


def relative_vectors(S, G):
    return S - G


def elevation_angle(G, GS):
    norm_G = np.linalg.norm(G, axis=1)
    norm_GS = np.linalg.norm(GS, axis=1)
    cos_phi = np.sum(G * GS, axis=1) / (norm_G * norm_GS)
    return np.rad2deg(np.arccos(np.clip(cos_phi, -1.0, 1.0)))


def horizontal_angle(G, GS):
    lat = np.deg2rad(lat_gs_deg)
    lon = np.deg2rad(lon_gs_deg)
    north = np.array([-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat)])
    east = np.array([-np.sin(lon), np.cos(lon), 0.0])
    up = np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])
    proj = GS - np.outer(np.dot(GS, up), up)
    x = np.dot(proj, east)
    y = np.dot(proj, north)
    angles = np.rad2deg(np.arctan2(x, y))
    return np.unwrap(np.deg2rad(angles), discont=np.pi) * 180.0 / np.pi


def compute_ctw(d, threshold):
    in_range = d <= threshold
    edges = np.diff(in_range.astype(int))
    starts = np.where(edges == 1)[0] + 1
    ends = np.where(edges == -1)[0] + 1
    if in_range[0]:
        starts = np.concatenate(([0], starts))
    if in_range[-1]:
        ends = np.concatenate((ends, [len(in_range)]))
    windows = [(s, e, (e - s) * dt) for s, e in zip(starts, ends)]
    return windows


def figure1_orbit_and_ground_station(S, G):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(S[:, 0], S[:, 1], S[:, 2], label='卫星轨道', color='tab:blue')
    ax.plot(G[:, 0], G[:, 1], G[:, 2], label='地面站轨迹', color='tab:orange')
    ax.scatter([0], [0], [0], color='k', s=20)
    ax.set_title('图1: 卫星轨道与地面站随地球自转运动')
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    ax.legend()
    ax.set_box_aspect([1, 1, 1])
    plt.tight_layout()
    return fig


def figure2_vector_geometry(S, G):
    idx = k_steps // 20
    Sk = S[idx]
    Gk = G[idx]
    GS_vec = Sk - Gk
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.quiver(0, 0, 0, Gk[0], Gk[1], Gk[2], color='tab:orange', length=np.linalg.norm(Gk), normalize=True)
    ax.quiver(Gk[0], Gk[1], Gk[2], GS_vec[0], GS_vec[1], GS_vec[2], color='tab:blue', length=np.linalg.norm(GS_vec), normalize=True)
    ax.scatter(*Gk, color='tab:orange', label='地面站位置')
    ax.scatter(*Sk, color='tab:blue', label='卫星位置')
    ax.set_title('图2: 卫星、地面站及参考矢量')
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    ax.legend()
    ax.set_box_aspect([1, 1, 1])
    plt.tight_layout()
    return fig


def figure3_distance_vs_time(d):
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(t / 3600.0, d, label='卫星-地面站距离')
    ax.axhline(d_max, color='tab:red', linestyle='--', label=f'最大通信范围 {d_max:.0f} km')
    ax.axvline(24, color='gray', linestyle=':')
    ax.axvline(48, color='gray', linestyle=':')
    ax.set_xlabel('时间 (h)')
    ax.set_ylabel('距离 (km)')
    ax.set_title('图3: 距离随时间变化')
    ax.legend()
    ax.grid(True)
    plt.tight_layout()
    return fig


def figure4_ground_station_perspective(S, G, d):
    lat = np.deg2rad(lat_gs_deg)
    lon = np.deg2rad(lon_gs_deg)
    north = np.array([-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat)])
    east = np.array([-np.sin(lon), np.cos(lon), 0.0])
    up = G / np.linalg.norm(G, axis=1)[:, None]
    GS = S - G
    x = np.dot(GS, east)
    y = np.dot(GS, north)
    fig, ax = plt.subplots(figsize=(8, 8))
    sc = ax.scatter(x, y, c=d, cmap='viridis', s=4)
    ax.set_title('图4: 地面站视角下卫星相对位置')
    ax.set_xlabel('East (km)')
    ax.set_ylabel('North (km)')
    ax.grid(True)
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label('距离 (km)')
    plt.tight_layout()
    return fig


def figure5_satellite_perspective(S, G, d):
    GS = G - S
    x_sat = GS[:, 0]
    y_sat = GS[:, 1]
    fig, axs = plt.subplots(1, 2, figsize=(14, 5))
    pass_mask = d <= d_max
    first_pass = np.where(pass_mask)[0]
    if len(first_pass) > 0:
        start = first_pass[0]
        end = first_pass[0] + 1
        while end < len(pass_mask) and pass_mask[end]:
            end += 1
        axs[0].scatter(x_sat[start:end], y_sat[start:end], c=d[start:end], cmap='plasma', s=10)
        axs[0].set_title('图5(a): 卫星视角单次通信通过')
    else:
        axs[0].set_title('图5(a): 没有检测到单次通过')
    axs[0].set_xlabel('Satellite frame X (km)')
    axs[0].set_ylabel('Satellite frame Y (km)')
    axs[0].grid(True)
    sc = axs[1].scatter(x_sat, y_sat, c=np.arange(len(t)), cmap='cool', s=2)
    axs[1].set_title('图5(b): 卫星视角三天距离向量轨迹')
    axs[1].set_xlabel('Satellite frame X (km)')
    axs[1].set_ylabel('Satellite frame Y (km)')
    axs[1].grid(True)
    fig.colorbar(sc, ax=axs[1], label='时间步')
    plt.tight_layout()
    return fig


def figure6_reference_vectors(S, G):
    idx = k_steps // 20
    Sk = S[idx]
    Gk = G[idx]
    GS = Sk - Gk
    lat = np.deg2rad(lat_gs_deg)
    lon = np.deg2rad(lon_gs_deg)
    north = np.array([-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat)])
    east = np.array([-np.sin(lon), np.cos(lon), 0.0])
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.quiver(*Gk, *north, length=1000, color='tab:green', label='北向参考')
    ax.quiver(*Gk, *GS, length=np.linalg.norm(GS), color='tab:blue', label='地面站->卫星')
    ax.quiver(*Gk, *east, length=1000, color='tab:red', label='东向参考')
    ax.scatter(*Gk, color='tab:orange', s=30, label='地面站')
    ax.scatter(*Sk, color='tab:purple', s=30, label='卫星')
    ax.set_title('图6: 地面站参考矢量与卫星方向')
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    ax.legend()
    ax.set_box_aspect([1, 1, 1])
    plt.tight_layout()
    return fig


def figure7_vertical_angle(phi):
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(t / 3600.0, phi, color='tab:purple')
    ax.set_title('图7: 垂直角随时间变化')
    ax.set_xlabel('时间 (h)')
    ax.set_ylabel('垂直角 φ (deg)')
    ax.grid(True)
    plt.tight_layout()
    return fig


def figure8_horizontal_angle(theta):
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(t / 3600.0, theta, color='tab:cyan')
    ax.set_title('图8: 水平角随时间变化')
    ax.set_xlabel('时间 (h)')
    ax.set_ylabel('水平角 θ (deg)')
    ax.grid(True)
    plt.tight_layout()
    return fig


def figure9_comm_range_geometry():
    fig, ax = plt.subplots(figsize=(8, 8))
    earth = plt.Circle((0, 0), RE, color='tab:blue', alpha=0.15, label='地球')
    orbit = plt.Circle((0, 0), RE + h, color='tab:gray', fill=False, linestyle='--', label='卫星轨道')
    ax.add_patch(earth)
    ax.add_patch(orbit)
    x_line = np.array([0.0, d_max])
    y_line = np.array([RE, 0.0])
    ax.plot([0, d_max], [RE, 0.0], color='tab:red', linestyle='-', label='最大通信切线')
    ax.scatter([0], [RE], color='tab:orange', label='地面站')
    ax.scatter([0], [RE + h], color='tab:purple', label='卫星在顶点')
    ax.set_aspect('equal', 'box')
    ax.set_xlim(-d_max * 1.1, d_max * 1.1)
    ax.set_ylim(RE - 1000, RE + h + 500)
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_title('图9: 最大通信范围几何图示')
    ax.legend()
    plt.tight_layout()
    return fig


def figure10_ctw(d):
    windows = compute_ctw(d, d_max)
    fig, axs = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    axs[0].plot(t / 3600.0, d, label='距离')
    for s, e, duration in windows:
        axs[0].axvspan(t[s] / 3600.0, t[e - 1] / 3600.0, color='tab:green', alpha=0.2)
    axs[0].axhline(d_max, color='tab:red', linestyle='--', label='最大通信范围')
    axs[0].set_title('图10(a): 三天 CTW 图形表示')
    axs[0].set_ylabel('距离 (km)')
    axs[0].legend()
    axs[0].grid(True)

    if windows:
        s, e, duration = windows[0]
        axs[1].plot(t[s:e] / 3600.0, d[s:e], color='tab:blue')
        axs[1].set_title('图10(b): 第一个 CTW')
        axs[1].set_xlabel('时间 (h)')
        axs[1].set_ylabel('距离 (km)')
        axs[1].grid(True)
    else:
        axs[1].text(0.5, 0.5, '未检测到 CTW', ha='center', va='center')
        axs[1].set_title('图10(b): 第一个 CTW')
        axs[1].set_xlabel('时间 (h)')
        axs[1].set_ylabel('距离 (km)')

    plt.tight_layout()
    return fig


def main():
    G = ground_station_positions()
    S = satellite_positions()
    GS = relative_vectors(S, G)
    d = np.linalg.norm(GS, axis=1)
    phi = elevation_angle(G, GS)
    theta = horizontal_angle(G, GS)

    figs = []
    figs.append(figure1_orbit_and_ground_station(S, G))
    figs.append(figure2_vector_geometry(S, G))
    figs.append(figure3_distance_vs_time(d))
    figs.append(figure4_ground_station_perspective(S, G, d))
    figs.append(figure5_satellite_perspective(S, G, d))
    figs.append(figure6_reference_vectors(S, G))
    figs.append(figure7_vertical_angle(phi))
    figs.append(figure8_horizontal_angle(theta))
    figs.append(figure9_comm_range_geometry())
    figs.append(figure10_ctw(d))

    output_dir = os.path.join(os.getcwd(), 'figures')
    os.makedirs(output_dir, exist_ok=True)

    for idx, fig in enumerate(figs, start=1):
        #fig.suptitle(f'生成图 {idx}')
        fig_path = os.path.join(output_dir, f'figure_{idx}.png')
        fig.savefig(fig_path, dpi=300, bbox_inches='tight')
        plt.close(fig)

    print(f'已将所有图像保存到: {output_dir}')


if __name__ == '__main__':
    main()
