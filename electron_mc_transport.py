import math
import numpy as np

# Constants

DEFAULT_NUMBER_DENSITY_M3 = 2.5e26     # default gas number density [particles / m^3]
EV_TO_J = 1.602176634e-19              # joule per eV
M_ELECTRON = 9.1093837015e-31          # kg


# Cross-section storage and loading


class XSecProcess:

    def __init__(
        self,
        kind,
        energy_v,
        xsec_v,
        loss_eV=0.0,
        name="",
        target_name="",
        mass_ratio=0.0,
    ):
        self.kind = kind.lower()
        self.energy_v = np.array(energy_v, dtype=float)
        self.xsec_v = np.array(xsec_v, dtype=float)
        self.loss_eV = float(loss_eV)
        self.name = name
        self.target_name = target_name
        self.mass_ratio = float(mass_ratio)

        order = np.argsort(self.energy_v)
        self.energy_v = self.energy_v[order]
        self.xsec_v = self.xsec_v[order]

    def value(self, energy_eV):
        """Return the cross section at energy_eV using linear interpolation."""
        if len(self.energy_v) == 0:
            return 0.0

        # No extrapolation outside the cross-section table.
        if energy_eV < self.energy_v[0] or energy_eV > self.energy_v[-1]:
            return 0.0

        return float(np.interp(energy_eV, self.energy_v, self.xsec_v))

    def __repr__(self):
        return "%s | target=%s | kind=%s | loss=%.3f eV | points=%d" % (
            self.name,
            self.target_name,
            self.kind,
            self.loss_eV,
            len(self.energy_v),
        )


def clean_target_name(name):
    """
    Try to get the target species name from a process name.

    Examples
    --------
    "He" -> "He"
    "He -> He(Singlet)" -> "He"
    "N2 -> N2*" -> "N2"
    """
    text = str(name).strip()

    if "<->" in text:
        text = text.split("<->")[0].strip()

    if "->" in text:
        text = text.split("->")[0].strip()

    words = text.split()
    if len(words) == 0:
        return ""

    return words[0]


def load_lxcat_table(filename, target_name=None, include_effective=False):
    """
    Load cross-section blocks from one LXCat-style file.

    The normal transport processes are:
        ELASTIC
        EXCITATION
        IONIZATION
        ATTACHMENT

    EFFECTIVE is ignored by default because it can double count processes
    if elastic, excitation, and ionization are also present.
    """
    processes = []

    valid_kinds = ["elastic", "excitation", "ionization", "attachment"]

    if include_effective:
        valid_kinds.append("effective")

    with open(filename, "r") as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        kind = line.lower()

        if kind in valid_kinds:
            name = lines[i + 1].strip()

            loss_eV = 0.0
            mass_ratio = 0.0

            if kind in ["elastic", "effective"]:
                try:
                    mass_ratio = float(lines[i + 2].split()[0])
                except (ValueError, IndexError):
                    mass_ratio = 0.0

            elif kind in ["excitation", "ionization"]:
                try:
                    loss_eV = float(lines[i + 2].split()[0])
                except (ValueError, IndexError):
                    loss_eV = 0.0

            if target_name is None:
                this_target_name = clean_target_name(name)
            else:
                this_target_name = target_name

            # Find the dashed line that starts this block's table.
            while i < len(lines) and not lines[i].strip().startswith("---"):
                i += 1

            # Read numbers until the dashed line that ends this block's table.
            i += 1
            energy_v = []
            xsec_v = []

            while i < len(lines) and not lines[i].strip().startswith("---"):
                words = lines[i].split()

                if len(words) >= 2:
                    try:
                        energy_v.append(float(words[0]))
                        xsec_v.append(float(words[1]))
                    except ValueError:
                        pass

                i += 1

            if len(energy_v) > 0:
                processes.append(
                    XSecProcess(
                        kind,
                        energy_v,
                        xsec_v,
                        loss_eV=loss_eV,
                        name=name,
                        target_name=this_target_name,
                        mass_ratio=mass_ratio,
                    )
                )

        i += 1

    if len(processes) == 0:
        raise RuntimeError("No usable cross-section blocks found in %s" % filename)

    return processes


def load_cross_section_table(filename, target_name=None):
    """Shorter name for load_lxcat_table()."""
    return load_lxcat_table(filename, target_name=target_name)


# ============================================================
# Material classes
# ============================================================

class TargetSpecies:

    def __init__(self, name, number_density_m3, processes):
        self.name = str(name)
        self.number_density_m3 = float(number_density_m3)
        self.processes = list(processes)

    def total_xsec(self, energy_eV):
        """Total microscopic cross section for this species."""
        total = 0.0

        for process in self.processes:
            total += process.value(energy_eV)

        return total

    def macro_xsec(self, energy_eV):
        """Macroscopic cross section for this species: Sigma = n sigma."""
        return self.number_density_m3 * self.total_xsec(energy_eV)

    def __repr__(self):
        return "%s | density=%.3e m^-3 | processes=%d" % (
            self.name,
            self.number_density_m3,
            len(self.processes),
        )


class TargetMaterial:
    """
    A target material made from one or more target species.

    For a pure gas, use one species.
    For a gas mixture, use several species with their own number densities.
    """

    def __init__(self, name, species):
        self.name = str(name)
        self.species = list(species)

        if len(self.species) == 0:
            raise ValueError("TargetMaterial needs at least one species.")

    def total_number_density(self):
        """Return total number density summed over all species."""
        total = 0.0

        for species in self.species:
            total += species.number_density_m3

        return total

    def total_macro_xsec(self, energy_eV):
        """Return total macroscopic cross section for the whole material."""
        total = 0.0

        for species in self.species:
            total += species.macro_xsec(energy_eV)

        return total

    def mean_free_path(self, energy_eV):
        """Mean free path lambda = 1 / Sigma_total."""
        sigma_macro_total = self.total_macro_xsec(energy_eV)

        if sigma_macro_total <= 0.0:
            return math.inf

        return 1.0 / sigma_macro_total

    def choose_interaction(self, energy_eV, null_factor=2.0):
        """
        Choose the target species and process for a real collision.

        The null_factor is used by the null-collision method.  With
        null_factor = 2, the total real-collision probability is about 1/2.
        """
        sigma_macro_total = self.total_macro_xsec(energy_eV)

        if sigma_macro_total <= 0.0:
            return None, None

        R = np.random.random()
        running = 0.0

        for species in self.species:
            for process in species.processes:
                sigma_macro = species.number_density_m3 * process.value(energy_eV)
                probability = sigma_macro / (null_factor * sigma_macro_total)
                running += probability

                if R < running:
                    return species, process

        return None, None

    def __repr__(self):
        return "%s | species=%d | density=%.3e m^-3" % (
            self.name,
            len(self.species),
            self.total_number_density(),
        )


def make_material_from_processes(
    material_name,
    target_name,
    number_density_m3,
    processes,
):
    """Make a one-species material from an already loaded process list."""
    species = TargetSpecies(target_name, number_density_m3, processes)
    return TargetMaterial(material_name, [species])


def load_lxcat_material(
    filename,
    material_name=None,
    target_name=None,
    number_density_m3=DEFAULT_NUMBER_DENSITY_M3,
):
    """
    Load one LXCat-style file and return a one-species TargetMaterial.
    """
    processes = load_lxcat_table(filename, target_name=target_name)

    if target_name is None:
        target_name = processes[0].target_name

    if material_name is None:
        material_name = target_name

    return make_material_from_processes(
        material_name,
        target_name,
        number_density_m3,
        processes,
    )


def make_gas_mixture(material_name, species_data):
    """
    Build a gas mixture.

    species_data should be a list of entries like:
        [name, number_density_m3, processes]

    Example
    -------
    air = make_gas_mixture("air", [
        ["N2", n_N2, n2_processes],
        ["O2", n_O2, o2_processes],
    ])
    """
    species_list = []

    for row in species_data:
        name = row[0]
        number_density_m3 = row[1]
        processes = row[2]

        species_list.append(TargetSpecies(name, number_density_m3, processes))

    return TargetMaterial(material_name, species_list)


# ============================================================
# Random directions
# ============================================================

def random_unit_vector():
    """Return a random 3D unit vector."""
    direction = np.random.normal(0.0, 1.0, 3)
    mag = np.linalg.norm(direction)

    if mag == 0.0:
        return np.array([1.0, 0.0, 0.0])

    return direction / mag


# ============================================================
# Electron class
# ============================================================

class ElectronParticle:
    """
    One electron moving through a target material.

    The update() method performs one Monte Carlo step:
        1. sample a collision time
        2. move the electron
        3. choose elastic / excitation / ionization / attachment / null
        4. update energy and direction
    """

    def __init__(
        self,
        energy_eV,
        material,
        x=0.0,
        y=0.0,
        z=0.0,
        t=0.0,
        idx=0,
        direction=None,
        energy_cutoff_eV=0.01,
        start_event="start",
    ):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.t = float(t)
        self.idx = idx

        self.energy_eV = float(energy_eV)
        self.material = material
        self.energy_cutoff_eV = float(energy_cutoff_eV)
        self.alive = True

        if direction is None:
            self.direction = random_unit_vector()
        else:
            self.direction = np.array(direction, dtype=float)
            mag = np.linalg.norm(self.direction)
            if mag == 0.0:
                self.direction = random_unit_vector()
            else:
                self.direction = self.direction / mag

        self.vx = 0.0
        self.vy = 0.0
        self.vz = 0.0

        self.distance_traveled = 0.0
        self.real_collision_count = 0
        self.history = []

        self.update_velocity_from_energy()
        self.record(start_event)

    def speed(self):
        """Return scalar speed."""
        return math.sqrt(self.vx * self.vx + self.vy * self.vy + self.vz * self.vz)

    def setenergy(self, energy_eV):
        """Set kinetic energy and update velocity magnitude."""
        self.energy_eV = max(float(energy_eV), 0.0)

        if self.energy_eV <= self.energy_cutoff_eV:
            self.alive = False

        self.update_velocity_from_energy()

    def update_velocity_from_energy(self):
        """Convert kinetic energy to speed using K = (1/2) m v^2."""
        if self.energy_eV <= self.energy_cutoff_eV:
            speed = 0.0
        else:
            energy_J = self.energy_eV * EV_TO_J
            speed = math.sqrt(2.0 * energy_J / M_ELECTRON)

        self.vx = speed * self.direction[0]
        self.vy = speed * self.direction[1]
        self.vz = speed * self.direction[2]

    def randomize_direction(self):
        """Scatter into a new random direction while keeping current energy."""
        self.direction = random_unit_vector()
        self.update_velocity_from_energy()

    def update_position(self, dt):
        """Move forward by dt seconds."""
        speed_before_move = self.speed()

        self.x += self.vx * dt
        self.y += self.vy * dt
        self.z += self.vz * dt
        self.t += dt

        self.distance_traveled += speed_before_move * dt

    def total_macro_xsec(self):
        """Total macroscopic cross section at the current electron energy."""
        return self.material.total_macro_xsec(self.energy_eV)

    def mean_free_path(self):
        """Mean free path lambda = 1 / Sigma_total."""
        return self.material.mean_free_path(self.energy_eV)

    def sample_dt(self):
        """
        Sample a null-collision Monte Carlo time step.

        This uses K = 2 * Sigma_total * v, so the null collision probability
        is about 1/2 when real cross sections are summed into Sigma_total.
        """
        sigma_macro_total = self.total_macro_xsec()
        v = self.speed()

        if sigma_macro_total <= 0.0 or v <= 0.0:
            self.alive = False
            return 0.0

        K = 2.0 * sigma_macro_total * v
        tau = 1.0 / K

        return np.random.exponential(tau)

    def choose_interaction(self):
        """
        Choose elastic, excitation, ionization, attachment, or null.
        """
        if self.speed() <= 0.0:
            return None, None

        return self.material.choose_interaction(self.energy_eV, null_factor=2.0)

    def scatter(self):
        """Elastic scatter: change direction, keep energy."""
        self.randomize_direction()
        self.real_collision_count += 1

    def excite(self, process):
        """Excitation: subtract excitation energy and scatter direction."""
        if self.energy_eV <= process.loss_eV:
            return False

        self.setenergy(self.energy_eV - process.loss_eV)
        self.randomize_direction()
        self.real_collision_count += 1

        return True

    def ionize(self, process, next_idx):
        """
        Ionization: create a secondary electron.

        The ionization threshold is removed first. The leftover energy is split
        randomly between the primary and secondary electrons.
        """
        if self.energy_eV <= process.loss_eV:
            return None, False

        leftover_energy = self.energy_eV - process.loss_eV
        fraction = np.random.random()

        primary_energy = fraction * leftover_energy
        secondary_energy = (1.0 - fraction) * leftover_energy

        self.setenergy(primary_energy)
        self.randomize_direction()
        self.real_collision_count += 1

        if secondary_energy <= self.energy_cutoff_eV:
            return None, True

        secondary = ElectronParticle(
            secondary_energy,
            self.material,
            x=self.x,
            y=self.y,
            z=self.z,
            t=self.t,
            idx=next_idx,
            direction=random_unit_vector(),
            energy_cutoff_eV=self.energy_cutoff_eV,
            start_event="birth",
        )

        return secondary, True

    def attach(self):
        """
        Attachment: remove the electron from the active simulation.
        """
        self.energy_eV = 0.0
        self.alive = False
        self.update_velocity_from_energy()
        self.real_collision_count += 1

    def update(self, next_idx=0):
        """
        Advance by one Monte Carlo update.

        Returns
        -------
        secondary : ElectronParticle or None
            New secondary electron from ionization, if created.
        event : str
            "scatter", "excitation", "ionization", "attachment", "null", or "dead".
        """
        if not self.alive:
            return None, "dead"

        dt = self.sample_dt()

        if not self.alive:
            self.record("dead")
            return None, "dead"

        self.update_position(dt)
        species, process = self.choose_interaction()
        secondary = None

        if process is None:
            event = "null"

        elif process.kind == "elastic":
            self.scatter()
            event = "scatter"

        elif process.kind == "excitation":
            happened = self.excite(process)
            event = "excitation" if happened else "null"

        elif process.kind == "ionization":
            secondary, happened = self.ionize(process, next_idx)
            event = "ionization" if happened else "null"

        elif process.kind == "attachment":
            self.attach()
            event = "attachment"

        else:
            event = "null"

        self.record(event)
        return secondary, event

    def record(self, event):
        """Store current state for printing and plotting."""
        self.history.append([self.t, self.x, self.y, self.z, self.energy_eV, event])

    def __repr__(self):
        return (
            "electron %d in %s: pos [%.3e, %.3e, %.3e] m, "
            "energy %.3f eV, speed %.3e m/s"
            % (
                self.idx,
                self.material.name,
                self.x,
                self.y,
                self.z,
                self.energy_eV,
                self.speed(),
            )
        )


# Backward-compatible names

ElectronHeParticle = ElectronParticle
HE_DENSITY_M3 = DEFAULT_NUMBER_DENSITY_M3
