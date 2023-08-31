"""Returns structures for tidy 3d simulation from gdsfactory components"""

 # Most of this code copied and stripped code from from gdsfactory\gplugins https://github.com/gdsfactory/gplugins.git 
from __future__ import annotations

import warnings

import gdsfactory as gf
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pydantic
import tidy3d as td
from gdsfactory.component import Component
from gdsfactory.components.extension import move_polar_rad_copy
from gdsfactory.config import logger
from gdsfactory.pdk import get_layer_stack
from gdsfactory.routing.sort_ports import sort_ports_x, sort_ports_y
from gdsfactory.technology import LayerStack
from gdsfactory.typings import ComponentSpec, Float2, Tuple
from tidy3d.plugins.mode import ModeSolver
from typing import List

from gplugins.tidy3d.materials import (
    get_index,
    get_medium,
)

def get_structures(
    component: ComponentSpec,
    layer_stack: LayerStack | None = None,
    clad_material: str = "sio2",
    center_wavelength: float = 1.55,
    material_name_to_tidy3d: None | dict[str, str] = None,
    is_3d: bool = True,
    **kwargs,
) -> List[td.Structure]:
    
    r"""Returns structures for tidy 3d simulation from gdsfactory components
    Args:
        component: gdsfactory Component.
        layer_stack: contains layer to thickness, zmin and material.
            Defaults to active pdk.layer_stack.
        clad_material: material for cladding.
        center_wavelength: center_wavelength
        material_name_to_tidy3d: dispersive materials have a wavelength
            dependent index. Maps layer_stack names with tidy3d material database names.
        is_3d: if False, does not consider Z dimension for faster simulations.
    keyword Args:
        medium: Background medium of simulation, defaults to vacuum if not specified.
        """
    
    component = gf.get_component(component)
    if not isinstance(component, Component):
        raise ValueError(f"component should be a gdsfactory.Component not {component}")
    
    layer_stack = layer_stack or get_layer_stack()
    layer_to_thickness = layer_stack.get_layer_to_thickness()
    layer_to_material = layer_stack.get_layer_to_material()
    layer_to_zmin = layer_stack.get_layer_to_zmin()
    layer_to_sidewall_angle = layer_stack.get_layer_to_sidewall_angle()

    material_name_to_tidy3d = material_name_to_tidy3d or {}

    if material_name_to_tidy3d:
        clad_material_spec = material_name_to_tidy3d[clad_material]
    else:
        clad_material_spec = (
            clad_material(center_wavelength) if callable(clad_material) else clad_material
        )

    clad = td.Structure(
        geometry=td.Box(
            size=(td.inf, td.inf, td.inf),
            center=(0, 0, 0),
        ),
        medium=get_medium(spec=clad_material_spec),
        name = 'background'
    )
    structures = [clad]

    if len(layer_to_thickness) < 1:
        raise ValueError(f"{component.get_layers()} not in {layer_to_thickness.keys()}")
    
    t_core = max(layer_to_thickness.values())
    component_layers = component.get_layers()

    
    for layer, thickness in layer_to_thickness.items():
        if layer in layer_to_material and layer in component_layers:
            zmin = layer_to_zmin[layer] if is_3d else -td.inf
            zmax = zmin + thickness if is_3d else td.inf

            material_name = layer_to_material[layer]

            if material_name in material_name_to_tidy3d:
                spec = material_name_to_tidy3d[material_name]
                medium = get_medium(spec=spec)
                index = get_index(spec=spec)
                logger.debug(
                    f"Add {layer}, {spec!r}, index = {index:.3f}, "
                    f"thickness = {thickness}, zmin = {zmin}, zmax = {zmax}"
                )
            else:
                material_index = (
                    material_name(center_wavelength)
                    if callable(material_name)
                    else material_name
                )
                medium = get_medium(material_index)

            polygons = td.PolySlab.from_gds(
                gds_cell=component._cell,
                gds_layer=layer[0],
                gds_dtype=layer[1],
                axis=2,
                slab_bounds=(zmin, zmax),
            )

            for polygon in polygons:
                geometry = td.Structure(geometry=polygon, medium=medium)
                structures.append(geometry)
    
    return structures

def get_monitors(
        component: ComponentSpec,
        port_margin: float = 0.5,
        monitor_offset: float = 0, 
        size_z: float = 3.,
        wavelength_start: float = 1.50,
        wavelength_stop: float = 1.60,
        wavelength_points: int = 10,
        is_3d = True,
        domain_size= [100,50,0],
        
)-> dict():
    cardinal_angles = [0., 90., 180., 270., 360.]
    
    # Add port monitors
    monitors = {}
    ports = sort_ports_x(sort_ports_y(component.get_ports_list()))
    for port in ports:
        port_name = port.name
        angle = port.orientation % 360

        # Find the cardinal angle with the smallest absolute difference
        nearest_angle = min(cardinal_angles, key=lambda cardinal: abs((angle - cardinal)))
        if nearest_angle % 180 == 0:
            theta = (angle - 180) % 360
        else:
            theta = (270 - angle) % 360

        width = port.width + 2 * port_margin
        size_x = width * abs(np.sin(nearest_angle * np.pi / 180))
        size_y = width * abs(np.cos(nearest_angle * np.pi / 180))
        size_x = 0 if size_x < 0.001 else size_x
        size_y = 0 if size_y < 0.001 else size_y
        size = (size_x, size_y, size_z)

        center = np.array(port.center).tolist() +[0]
        center_offset = move_polar_rad_copy(
            np.array(port.center), angle=angle * np.pi / 180, length= monitor_offset
        )
        center_offset = center_offset.tolist() + [0] 
        

        wavelengths = np.linspace(wavelength_start, wavelength_stop, wavelength_points)
        freqs = td.constants.C_0 / wavelengths
        freq0 = td.constants.C_0 / np.mean(wavelengths)
        fwidth = freq0 / 10

        monitors[port_name] = td.ModeMonitor(
            center=center_offset,
            size=size,
            freqs=freqs,
            mode_spec=td.ModeSpec(num_modes=2,angle_theta=np.deg2rad(theta)),
            name=f"mode_{port.name}",
        )

    zcenter = 0.1 if is_3d else 0
    domain_monitor = td.FieldMonitor(
        center=[0, 0, zcenter],
        size=domain_size if is_3d else [td.inf, td.inf, 0],
        freqs=[freq0],
        name="field",
    )
    monitors = list(monitors.values())
    monitors += [domain_monitor]
    return monitors

def get_sources(
        component: ComponentSpec,
        port_source_names: list(str) = ["o1"],
        source_offset: float = 0.,
        port_margin: float = 0.5,
        wavelength_start: float = 1.50,
        wavelength_stop: float = 1.60,
        wavelength_points: int = 10,
        size_z: float = 3.,
        is_3d: bool = True,
        num_modes: int = 1,
        mode_index: int = 0,

)-> dict():
    
    wavelengths = np.linspace(wavelength_start, wavelength_stop, wavelength_points)
    freqs = td.constants.C_0 / wavelengths
    freq0 = td.constants.C_0 / np.mean(wavelengths)
    fwidth = freq0 / 10
    #sources={}
    for port_source_name in port_source_names:    
        if port_source_name not in component.ports:
            warnings.warn(
                f"port_source_name={port_source_name!r} not in {list(component.ports.keys())}"
            )
            port_source = component.get_ports_list(port_type="optical")[0]
            port_source_name = port_source.name
            warnings.warn(f"Selecting port_source_name={port_source_name!r} instead.")

        # Add source
        port = component.ports[port_source_name]
        name=port.name
        angle= port.orientation % 360
        cardinal_angles = [0., 90., 180., 270., 360.]

        # Find the cardinal angle with the smallest absolute difference
        nearest_angle = min(cardinal_angles, key=lambda cardinal: abs((angle - cardinal)))
        if nearest_angle % 180 == 0:
            theta = (angle - 180) % 360
        else:
            theta = (270 - angle) % 360

        width = port.width + 2 * port_margin
        size_x = width * abs(np.sin(nearest_angle * np.pi / 180))
        size_y = width * abs(np.cos(nearest_angle * np.pi / 180))

        size_x = 0 if size_x < 0.001 else size_x
        size_y = 0 if size_y < 0.001 else size_y
        size_z = size_z if is_3d else td.inf

        source_size = [size_x, size_y, size_z]
        source_center = port.center.tolist() + [0]  # (x, y, z=0)

        xy_shifted = move_polar_rad_copy(
            np.array(port.center), angle=angle * np.pi / 180, length=source_offset
        )
        source_center_offset = xy_shifted.tolist() + [0]  # (x, y, z=0)

    #  theta=theta if int(angle) in {45., 225.} else -theta
        mode_spec=td.ModeSpec(num_modes= num_modes, angle_theta=np.deg2rad(theta))

        sources = td.ModeSource(
            size=source_size,
            center=source_center_offset,
            mode_spec=mode_spec,
            source_time=td.GaussianPulse(freq0=freq0, fwidth=fwidth),
            mode_index=mode_index,
            direction="+", 
            name = port_source_name,
            num_freqs= 5,
        )
        return sources

def parse_port_eigenmode_coeff(
    port_name: str, component: ComponentSpec, sim_data: td.SimulationData, mode_index: int=0
) -> Tuple[np.ndarray]:
    """Given a port and eigenmode coefficient result, returns the coefficients \
    relative to whether the wavevector is entering or exiting simulation.

    """
    sim = sim_data.simulation
    orientation = component.ports[port_name].orientation
    direction_inp = "-"
    direction_out = "-"
    coeff_inp = sim_data.monitor_data[f"mode_{port_name}"].amps.sel(direction = direction_inp, mode_index=mode_index)
    coeff_out = sim_data.monitor_data[f"mode_{port_name}"].amps.sel(direction = direction_out, mode_index=mode_index)
    return coeff_inp.values.flatten(), coeff_out.values.flatten()

def get_wavelengths(port_name: str, sim_data: td.SimulationData) -> np.ndarray:
    coeff_inp = sim_data.monitor_data[f"mode_{port_name}"].amps.sel(direction="+")
    freqs = coeff_inp.f
    return td.constants.C_0 / freqs.values

def get_sparam(port_out: str, component: ComponentSpec, sim_data: td.SimulationData, mode_index = 0
) -> np.ndarray:
    sim = sim_data.simulation
    port_name_source = sim.sources[0].name
    source_entering, source_exiting = parse_port_eigenmode_coeff(
            port_name=port_name_source, component=component, sim_data=sim_data, mode_index = mode_index
        )
    
    monitor_entering, monitor_exiting = parse_port_eigenmode_coeff(
                port_name=port_out,  component=component, sim_data=sim_data, mode_index = mode_index
            )
    sij = monitor_exiting / source_entering


    return sij
    
