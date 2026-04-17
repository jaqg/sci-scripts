import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# +---------------------------------------------------------------------------+
# |                                CLASSES                                    |
# +---------------------------------------------------------------------------+
class SpectrumData:
    """
    Clase para lectura y conversion de datos
    """
    H = 4.135667696e-15  # eV s
    C_M = 2.99792458e8    # m s-1
    C_CM = 2.99792458e10  # cm s-1
    C_NM = 2.99792458e17  # nm s-1

    def __init__(self, filename, rows_to_skip=0):
        self.filename = os.path.join(os.path.dirname(__file__), filename)
        self.energy, self.me = self.read_data(rows_to_skip)
        # self.wl = self.wavelength_to_energy(self.energy)
        self.cm_inv = self.energy_to_inv_cm()
    
    def read_data(self, rows_to_skip):
        return np.loadtxt(self.filename, unpack=True, skiprows=rows_to_skip)
    
    # def wavelength_to_energy(self):
    #     return self.H * self.C_M / self.wl
    
    def energy_to_inv_cm(self):
        # return self.energy / (self.H * self.C_CM)
        return self.energy * 8065.73
    
    def energy_to_nm(self):
        return self.H * self.C_NM / self.energy

    @staticmethod
    def load_from_folder(folder='data'):
        """
        Lee todos los archivos .dat de la carpeta especificada y devuelve un
        diccionario con los objetos de SpectrumData basandose en el nombre del archivo.
        Nota: los caracteres '-' se reemplazan por '_'
        """
        folder_path = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), folder)
        LSs = {}
        for file in os.listdir(folder_path):
            if file.endswith(".dat"):
                file_path = os.path.join(folder_path, file)
                var_name = os.path.splitext(file)[0].replace('-', '_')
                LSs[var_name] = SpectrumData(file_path)
        return LSs

class SpectrumPlotter:
    """
    Clase para generar graficos
    """
    def __init__(self, title, plot_all_TF=True):
        self.title = title
        self.plot_all_TF = plot_all_TF
        self.data_sets = []
        self.fig = None

    def set_figure(self, fig):
        self.fig = fig
    
    def add_data(self, spectrum_data, label, color=None, linestyle=None):
        self.data_sets.append((spectrum_data, label, color, linestyle))
    
    def plot(self, ax):
        for data, label, color, linestyle in self.data_sets:
            # ax.plot(data.cm_inv, data.me, marker='None', label=label, color=color, linestyle=linestyle)
            ax.plot(data.energy, data.me, marker='None', label=label, color=color, linestyle=linestyle)
        ax.legend()
    
    def configure_plot(self, ax, xlim_min=5, xlim_max=6, no_xlabel=False, no_ylabel=False):
        # ax.set_xlim(40000, 50000)  # cm-1
        ax.set_xlim(xlim_min, xlim_max)  # eV
        if no_xlabel:
            x_label=r''
        else:
            x_label=r'Energy (eV)'
            # xlabel=r'$\tilde{\nu}\ (\mathrm{cm^{-1}})$',
        if no_ylabel:
            y_label=r''
        else:
            # ylabel=r'$\Delta \varepsilon \left( \mathrm{L\,mol^{-1}\,cm^{-1}} \right)$'
            y_label=r'Lineshape (a.u.)'
        ax.set(
            title=self.title,
            xlabel=x_label,
            ylabel=y_label,
        )
    
    def save_figure(self, filename):
        if self.fig is None:
            raise RuntimeError("Figure not assigned. Use set_figure() before saving.")

        fn = os.path.dirname(os.path.abspath(sys.argv[0]))+'/{}'.format(filename)
        self.fig.savefig(fn, transparent=True, bbox_inches='tight')
        for ext in ['pdf']:
            fn_e = filename.replace('.pdf', f'.{ext}')
            fn = os.path.dirname(os.path.abspath(sys.argv[0]))+'/{}'.format(fn_e)
            self.fig.savefig(fn)
