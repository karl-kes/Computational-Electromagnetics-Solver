import pandas as pd
import numpy as np
import plotly.graph_objects as go
import glob
import os

def load_volume_data(file_path):
    cols = ['x', 'y', 'z', 'Fx', 'Fy', 'Fz', 'mag']
    df = pd.read_csv(file_path, names=cols)
    
    nx, ny, nz = int(df['x'].max() + 1), int(df['y'].max() + 1), int(df['z'].max() + 1)
    
    vol = df.sort_values(by=['z', 'y', 'x'])['mag'].values.reshape((nz, ny, nx))
    return vol, (nx, ny, nz)

# 1. Load data files
files = sorted(glob.glob("output/E*.csv"), 
               key=lambda x: int(''.join(filter(str.isdigit, x))))

volume_list = []
for f in files:
    v, shape = load_volume_data(f)
    volume_list.append(np.sqrt(v))

nz, ny, nx = shape
z, y, x = np.mgrid[0:nz, 0:ny, 0:nx]

# 2. Build Animation Frames with Dynamic Scaling
frames = []
for i, vol in enumerate(volume_list):
    v_max = vol.max() if vol.max() > 1e-15 else 1.0
    v_min = v_max * 0.05 
    
    frames.append(go.Frame(
        data=[go.Volume(
            value=vol.flatten(),
            isomin=v_min,
            isomax=v_max
        )],
        name=str(i)
    ))

# 3. Create Final Figure
first_max = volume_list[0].max() if volume_list[0].max() > 1e-15 else 1.0
fig = go.Figure(
    data=[go.Volume(
        x=x.flatten(), y=y.flatten(), z=z.flatten(),
        value=volume_list[0].flatten(),
        isomin=first_max * 0.05,
        isomax=first_max,
        opacity=0.15,
        surface_count=30,
        colorscale='Inferno',
        showscale=True,
        colorbar=dict(title="Field Intensity", thickness=20),
        caps=dict(x_show=False, y_show=False, z_show=False)
    )],
    layout=go.Layout(
        template="plotly_dark",
        title="Field Propogation",
        scene=dict(
            xaxis=dict(title="X"),
            yaxis=dict(title="Y"),
            zaxis=dict(title="Z"),
            aspectmode='cube'
        ),
        updatemenus=[dict(
            type="buttons",
            buttons=[
                dict(label="Play", method="animate", 
                    args=[None, {"frame": {"duration": 40, "redraw": True}, "fromcurrent": True}]),
                dict(label="Pause", method="animate", 
                     args=[[None], {"mode": "immediate"}])
            ],
            x=0.1, y=0
        )]
    ),
    frames=frames
)

fig.show()