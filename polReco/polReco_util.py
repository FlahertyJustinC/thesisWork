import os
import pickle
import numpy as np
import scipy.signal
from pyrex.internal_functions import (normalize, complex_bilinear_interp, complex_interp)
# import pyrex.custom.analysis as analysis
# import pyrex.custom.araroot as araroot
import pyrex.custom.envelope_reco as reco
import pyrex.custom.ara as ara
from pyrex.internal_functions import (normalize, complex_bilinear_interp,
                                      complex_interp)

"""
Some of these functions were modified from PyREx
"""


def load_tof_grid_data(filename, antenna_positions=None):
    with np.load(filename) as f:
        if antenna_positions is not None:
            for i, pos in enumerate(antenna_positions):
                if not np.all(np.isclose(pos, f['antennas'][i])):
                    raise ValueError("File does not match antenna positions")
        rs = f['rs']
        thetas = f['thetas']
        phis = f['phis']
        tofs_dir = f['tofs_dir']
        tofs_ref = f['tofs_ref']
    return (tofs_dir, tofs_ref), (rs, thetas, phis)

def calculate_interferometric_grid(waveforms, tofs, hilbert=False):
    from itertools import combinations
    bins = len(waveforms[0].values)
    dt = waveforms[0].dt
    t0s = np.asarray([wave.times[0] for wave in waveforms])
    xtof = tofs[:, :, :, :, np.newaxis] - tofs[:, :, :, np.newaxis, :]
    xtdiff = t0s[:, np.newaxis] - t0s[np.newaxis, :]
    shifts = (xtof-xtdiff)/dt
    int_shifts = np.asarray(shifts, dtype=np.int_) + bins - 1
    nan_shifts = (np.isnan(shifts) | np.asarray(int_shifts<0)
                  | np.asarray(int_shifts>=2*bins-1))
    int_shifts[nan_shifts] = 0
    xcorrs = reco.calculate_cross_correlations(waveforms, hilbert=hilbert)
    inter_data = np.zeros(tofs.shape[:3])
    combos = list(combinations(range(len(waveforms)), 2))
    for idx1, idx2 in combos:
        val = xcorrs[idx1, idx2, int_shifts[:, :, :, idx1, idx2]]
        val[nan_shifts[:, :, :, idx1, idx2]] = 0
        inter_data += val
    inter_data /= len(combos)


    return inter_data

def doReco(antenna_waveforms, plot_map=False, pol = 0):
    if(pol==0):
        tof_data, grid_points = load_tof_grid_data("tofs_ara02_vpols_3000m_spherical.npz")
    else:
        tof_data, grid_points = load_tof_grid_data("tofs_ara02_hpols_3000m_spherical.npz")

    # antenna_waveforms = pyrex_array
    mask = None
    if mask is not None:
        if len(mask)!=len(antennas):
            raise ValueError("Mask must have same length as antenna list")
        antenna_waveforms = [pos for pos, use in zip(antenna_waveforms, mask)
                             if use]
    else:
        mask = slice(None)
    if(pol==0):
        antenna_waveforms = antenna_waveforms[:8]
    else:
        antenna_waveforms = antenna_waveforms[8:]

    # Get the best interferometric fit among the time of flight data sets
    inter_max = 0
    best_data = None
    for tofs in tof_data:
        inter_data = reco.calculate_interferometric_grid(antenna_waveforms,
                                                    tofs[:, :, :, mask],
                                                    hilbert=True)
        max_val = np.max(inter_data)
        if best_data is None or max_val>inter_max:
            best_data = np.array(inter_data, copy=True)
            inter_max = max_val
    # Get the vertex of the highest interferometric value
    i, j, k = np.unravel_index(np.argmax(best_data), best_data.shape)
    inter_vtx = np.array((grid_points[0][i],
                          grid_points[1][j],
                          grid_points[2][k]))
    # print(inter_data[1][0])
    # print(inter_vtx, inter_max)
    # plot_map = False
    if(plot_map):
        max_idx = np.unravel_index(np.argmax(inter_data), inter_data.shape)
        # print("Hilbert Coherence", name+":", np.max(inter_data))
        with np.load(os.path.join('tofs_ara02_vpols_calpulser_spherical.npz')) as f:
            antenna_positions = None
            if antenna_positions is not None:
                for i, pos in enumerate(antenna_positions):
                    if not np.all(np.isclose(pos, f['antennas'][i])):
                        raise ValueError("File does not match antenna positions")
            rs = f['rs']
            thetas = f['thetas']
            phis = f['phis']
            tofs_dir = f['tofs_dir']
            tofs_ref = f['tofs_ref']
        tof_data = (tofs_dir, tofs_ref)
        grid_points = (rs, thetas, phis)
        mesh_thetas = np.concatenate((thetas, [thetas[-1]+(thetas[1]-thetas[0])]))
        mesh_phis = np.concatenate((phis, [phis[-1]+(phis[1]-phis[0])]))

        for name, tofs in zip(('Direct', 'Reflected'), tof_data):
            inter_data = reco.calculate_interferometric_grid(antenna_waveforms, tofs, hilbert=True)
            max_idx = np.unravel_index(np.argmax(inter_data), inter_data.shape)
            print("Hilbert Coherence", name+":", np.max(inter_data))
            plt.pcolormesh(mesh_phis, mesh_thetas, inter_data[0, :, :], vmin=0, vmax=1)
            plt.colorbar()
            plt.scatter(phis[max_idx[2]], thetas[max_idx[1]], color='k', marker=',', s=1)
            plt.title("Fullband Hilbert Reconstruction "+name)
            plt.show()

    return inter_vtx, inter_max

def interpolate_filter(frequencies):
        """
        Generate interpolated filter values for given frequencies.
        Calculate the interpolated values of the antenna system's filter gain
        data for some frequencies.
        Parameters
        ----------
        frequencies : array_like
            1D array of frequencies (Hz) at which to calculate gains.
        Returns
        -------
        array_like
            Complex filter gain in voltage for the given `frequencies`.
        """
        ARAfilter = ara.antenna.ALL_FILTERS_DATA
        filt_response = ARAfilter[0]
        filt_freqs = ARAfilter[1]
        return complex_interp(
            x=frequencies, xp=filt_freqs, fp=filt_response,
            method='euler', outer=0
        )

def directional_response(theta, phi, polarization=np.array([0,0,1])):
        """
        Generate the (complex) frequency-dependent directional response.
        For given angles and polarization direction, use the model of the
        directional and polarization gains of the antenna to generate a
        function for the interpolated response of the antenna with respect to
        frequency. Used with the `frequency_response` method to calculate
        effective heights.
        Parameters
        ----------
        theta : float
            Polar angle (radians) from which a signal is arriving.
        phi : float
            Azimuthal angle (radians) from which a signal is arriving.
        polarization : array_like
            Normalized polarization vector in the antenna coordinate system.
        Returns
        -------
        function
            A function which returns complex-valued voltage gains for given
            frequencies, using the values of incoming angle and polarization.
        See Also
        --------
        ARAAntenna.frequency_response : Calculate the (complex) frequency
                                        response of the antenna.
        """
#         e_theta = [np.cos(theta) * np.cos(phi),
#                    np.cos(theta) * np.sin(phi),
#                    -np.sin(theta)]
#         e_phi = [-np.sin(phi), np.cos(phi), 0]
#         theta_factor = np.dot(polarization, e_theta)
#         phi_factor = np.dot(polarization, e_phi)
        theta_factor = 1
        phi_factor = 1
        theta_gains = complex_bilinear_interp(
            x=np.degrees(theta), y=np.degrees(phi),
            xp=response_zens,
            yp=response_azis,
            fp=theta_response,
            method='cartesian'
        )
        phi_gains = complex_bilinear_interp(
            x=np.degrees(theta), y=np.degrees(phi),
            xp=response_zens,
            yp=response_azis,
            fp=phi_response,
            method='cartesian'
        )
        freq_interpolator = lambda frequencies: complex_interp(
            x=frequencies, xp=response_freqs,
            fp=theta_factor*theta_gains + phi_factor*phi_gains,
            method='euler', outer=0
        )
        return freq_interpolator

def frequency_response(frequencies):
    """
    Calculate the (complex) frequency response of the antenna.
    Rather than handling the entire frequency response of the antenna, this
    method is being used to convert the frequency-dependent gains from the
    `directional_response` method into effective heights.
    Parameters
    ----------
    frequencies : array_like
        1D array of frequencies (Hz) at which to calculate gains.
    Returns
    -------
    array_like
        Complex gains in voltage for the given `frequencies`.
    See Also
    --------
    ARAAntenna.directional_response : Generate the (complex) frequency
                                      dependent directional response.
    """
    # From AraSim GaintoHeight function, with gain calculation moved to
    # the directional_response method.
    # gain=4*pi*A_eff/lambda^2 and h_eff=2*sqrt(A_eff*Z_rx/Z_air)
    # Then 0.5 to calculate power with heff (cancels 2 above)
    heff = np.zeros(len(frequencies))
    # The index of refraction in this calculation should be the index of
    # the ice used in the production of the antenna model.
    n = 1.78
    heff[frequencies!=0] = np.sqrt((3e8/frequencies[frequencies!=0]/n)**2
                                   * n*50/377 /(4*np.pi))
    return heff

def doFFT(time, volts):
    """
    Calculate the Fast-Fourier transform (FFT) of a signal.
    ----------
    time : array_like
        1D array of times (ns).
    volts : array_like
        1D array of amplitudes (mV).
    Returns
    -------
    fft : array_like
        Amplitude in f-domain.
    freq : array_like
        Frequencies in MHz
    """
    fft = scipy.fft.rfft(np.array(volts))
    dT = abs(time[1]-time[0])
    freq = 1000*scipy.fft.rfftfreq(n=len(time), d=dT)
    return fft, freq, dT

def doInvFFT(spectrum):
    """
    Calculate the inverse Fast-Fourier transform (FFT) of a signal.
    ----------
    spectrum : array_like
        1D array of amplitudes in f-domain
    Returns
    -------
    fft_i_v : array_like
        Amplitudes in mV.
    """
    fft_i_v= scipy.fft.irfft(spectrum)
    return fft_i_v

def deDisperse_filter(time, voltage):
    """
    Apply inverse of ARA filter response phase (amplitudes remain the same)
    ----------
    time : array_like
        1D array of times (ns)
    voltage : array_like
        1D array of amplitudes (mV).

    Returns
    -------
    time : array_like
        1D array of times (ns)
    deDis_wf : array_like
        1D array of amplitudes (mV) of de-dispersed waveform.
    """
    fft_v, fft_f, dT = doFFT(time,voltage)
    response = np.array(interpolate_filter(fft_f*1E6))
    response = np.divide(response,abs(response))
    deDis_wf = np.divide(fft_v,response)
    deDis_wf = np.nan_to_num(deDis_wf)
    deDis_wf = doInvFFT(deDis_wf)
    return time, deDis_wf

def deDisperse_antenna(time, voltage, theta, phi):
    """
    Apply inverse of ARA antenna response phase (amplitudes remain the same)
    ----------
    time : array_like
        1D array of times (ns)
    voltage : array_like
        1D array of amplitudes (mV).
    theta : double
        Incoming signal theta direction
    phi : double
        Incoming signal phi direction

    Returns
    -------
    time : array_like
        1D array of times (ns)
    deDis_wf : array_like
        1D array of amplitudes (mV) of de-dispersed waveform.
    """
    fft_v, fft_f, dT = doFFT(time, voltage)
    dir_res = directional_response(theta,phi)(fft_f*1E6)
    heff = dir_res * frequency_response(fft_f*1E6)
    response = dir_res*heff
    response = np.divide(response,abs(response))
    deDis_wf = np.divide(fft_v,response)
    deDis_wf = np.nan_to_num(deDis_wf)
    deDis_wf = doInvFFT(deDis_wf)
    return time, deDis_wf

def deDisperse(time, voltage, theta, phi):
    """
    Apply inverse of ARA antenna+filter response phase (amplitudes remain the same)
    ----------
    time : array_like
        1D array of times (ns)
    voltage : array_like
        1D array of amplitudes (mV).
    theta : double
        Incoming signal theta direction
    phi : double
        Incoming signal phi direction

    Returns
    -------
    time : array_like
        1D array of times (ns)
    deDis_wf : array_like
        1D array of amplitudes (mV) of de-dispersed waveform.
    """
    sampRate = len(time)/(max(time)-min(time))
    b,a = signal.bessel(4, [0.15,0.4], 'bandpass', analog=False, fs=sampRate)
    voltage = signal.lfilter(b, a, voltage)
    fft_v, fft_f, dT = doFFT(time,voltage)
    response_filter = np.array(interpolate_filter(fft_f*1E6))
    dir_res = directional_response(theta,phi)(fft_f*1E6)
    heff = dir_res * frequency_response(fft_f*1E6)
    response_antenna = dir_res*heff
    response = response_filter + response_antenna
    response = np.divide(response,abs(response))
    deDis_wf = np.divide(fft_v,response)
    deDis_wf = np.nan_to_num(deDis_wf)
    deDis_wf = doInvFFT(deDis_wf)
    return time, deDis_wf

def deConvolve_antenna(time, voltage, theta, phi, pol_ant, channel=None, station=None, configuration=None):
    """
    Apply inverse of ARA antenna response
    ----------
    time : array_like
        1D array of times (ns)
    voltage : array_like
        1D array of amplitudes (mV).

    theta, phi, pol_ant: floats
    theta_antenna (radians), phi_antenna (radians), pol_antenna [0:vpol, 1:hpol]
    Returns
    -------
    time : array_like
        1D array of times (ns)
    deDis_wf : array_like
        1D array of amplitudes (mV) of de-convolved waveform.
    """
    import scipy.signal as signal
    polarization=np.array([-np.sin(phi),np.cos(phi),-1/np.sin(theta)]) #This is just a factor to cancel out aspects of the antenna response in PyRex. - JCF 8/30/2022
    if(pol_ant == 0):
        ant = ara.VpolAntenna(name="Dummy Vpol", position=(0, 0, 0), power_threshold=0, channel=channel, station=station, configuration=configuration)
        # ant.set_orientation(z_axis=(0, 0, 1), x_axis=(1, 0, 0))#Adding to convert from global coordinates to local antenna coords.
    elif(pol_ant == 1):
        ant = ara.HpolAntenna(name="Dummy Hpol", position=(0, 0, 0), power_threshold=0, channel=channel, station=station, configuration=configuration)
        # ant.set_orientation(z_axis=(0, 0, 1), x_axis=(1, 0, 0))

    sampRate = len(time)/(max(time)-min(time))
    b,a = signal.bessel(4, [0.15,0.4], 'bandpass', analog=False, fs=sampRate)
    fft_v, fft_f, dT = doFFT(time,voltage)
    # response_filter = np.array(interpolate_filter(fft_f*1E6))
    # response_filter = np.array(ara.interpolate_filter(fft_f*1E6))
    dir_res = ant.antenna.directional_response(theta=theta, phi=phi, polarization=polarization)(fft_f*1E6)
    heff = ant.antenna.frequency_response(fft_f*1E6)
    response_antenna = dir_res*heff
    response = response_antenna
    deDis_wf = np.divide(fft_v,abs(response))
    response = np.divide(response,abs(response))
    deDis_wf = np.divide(deDis_wf,response)  #What is the purpose of this step?  The waveform seems to break without it. - JCF 6/19/2023
    deDis_wf = np.nan_to_num(deDis_wf)
    revert = doInvFFT(deDis_wf)
    deDis_wf = signal.lfilter(b, a, revert)
    # deDis_wf = revert
    if (len(time) > len(deDis_wf)):
        time = time[:-1]
    return time, deDis_wf
    #vetted!


def PolAngleStokes(Hpol,Vpol):
    return np.degrees(0.5*np.arctan2(2*Hpol*Vpol,(Vpol**2-Hpol**2)))
def PolRatio(Hpol,Vpol):
    return np.degrees(np.arctan2(Hpol,Vpol))

def getRFChannel(string, channel): #mapping from AraSim to RF channel chain
    RFChannel = 0
    if(string == 0):
        if(channel == 0):
            RFChannel = 5
        elif(channel == 1):
            RFChannel = 13
        elif(channel == 2):
            RFChannel = 1
        elif(channel == 3):
            RFChannel = 9

    elif(string == 1):
        if(channel == 0):
            RFChannel = 6
        elif(channel == 1):
            RFChannel = 14
        elif(channel == 2):
            RFChannel = 2
        elif(channel == 3):
            RFChannel = 10

    elif(string == 2):
        if(channel == 0):
            RFChannel = 7
        elif(channel == 1):
            RFChannel = 15
        elif(channel == 2):
            RFChannel = 3
        elif(channel == 3):
            RFChannel = 11

    if(string == 3):
        if(channel == 0):
            RFChannel = 4
        elif(channel == 1):
            RFChannel = 12
        elif(channel == 2):
            RFChannel = 0
        elif(channel == 3):
            RFChannel = 8

    return int(RFChannel)

def PolVectorReco(Peak_V, Peak_H, theta, phi):
    R = (Peak_H/Peak_V)
    denom = np.sqrt(1+R**2)
    Px = -(np.cos(theta)*np.cos(phi)-R*np.sin(phi))/denom
    Py = -(R*np.cos(phi)+np.cos(theta)*np.sin(phi))/denom
    Pz = np.sin(theta)/denom
    np.set_printoptions(suppress=True)
    # if Peak_V>0:
    #     Px = -Px
    #     Py = -Py
    #     Pz = -Pz
    return np.array([Px,Py,Pz])

def findMaxSign(s1):
    if(abs(max(s1))>=abs(min(s1))):
        value = max(s1)
    else:
        value = min(s1)
    return value

def findAveragePeak(s1):
    posPeak = max(s1)
    negPeak = min(s1)
    return (posPeak-negPeak)/2

def getResponseAraSim(theta, phi, freq, pol):
    """
    Get antenna response from AraSim
    ----------
    theta : double
        Arriving angle at antenna (radians)
    phi : double
        Arriving angle at antenna (radians)
    freq : array_like
        1D array of frequencies (MHz)
    pol : integer
        0: Vpol, 1: Hpol

    Returns
    -------
    freq : array_like
        1D array of input frequencies (Hz) WHY???
    heffs : array_like
        1D array of complex effective heights
    """
    from ROOT import TCanvas, TGraph
    from ROOT import gROOT
    import ROOT
    import os
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    from ROOT import gInterpreter, gSystem
    from ROOT import TChain, TSelector, TTree
    import cmath

    gInterpreter.ProcessLine('#include "/cvmfs/ara.opensciencegrid.org/trunk/centos7/source/AraSim/Position.h"')
    gInterpreter.ProcessLine('#include "/cvmfs/ara.opensciencegrid.org/trunk/centos7/source/AraSim/Report.h"')
    gInterpreter.ProcessLine('#include "/cvmfs/ara.opensciencegrid.org/trunk/centos7/source/AraSim/Detector.h"')
    gInterpreter.ProcessLine('#include "/cvmfs/ara.opensciencegrid.org/trunk/centos7/source/AraSim/Settings.h"')

    gSystem.Load('/cvmfs/ara.opensciencegrid.org/trunk/centos7/source/AraSim/libAra.so') #load the simulation event library. You might get an error asking for the eventSim dictionry. To solve that, go to where you compiled AraSim, find that file, and copy it to where you set LD_LIBRARY_PATH.


    file_list = []
    file_list.append("/cvmfs/ara.opensciencegrid.org/trunk/centos7/source/AraSim/outputs/AraOut.default_A2_c1_E610_readIn.txt.runAraSim_comparison_input_test_1E19_2.txt.root")

    simSettingsTree = TChain("AraTree")
    simTree = TChain("AraTree2")

    for line in file_list:
        simTree.AddFile(line)
        simSettingsTree.AddFile(line)

    reportPtr = ROOT.Report()
    detectorPtr = ROOT.Detector()

    simTree.SetBranchAddress("report", ROOT.AddressOf(reportPtr))
    simSettingsTree.SetBranchAddress("detector", ROOT.AddressOf(detectorPtr))
    numEvents = simTree.GetEntries()


    simTree.GetEntry(0)
    simSettingsTree.GetEntry(0)

    theta = np.degrees(theta)
    phi = np.degrees(phi)
    freq = freq*1E6

    # dt = 0.3125e-9 # seconds
    # ff = np.fft.rfftfreq(int(1280/2), dt)
    gains = []
    heffs = []
    filter_gains = []
    # phases = []
    for f in freq:
        gain = detectorPtr.GetGain_1D_OutZero(f/1e6, theta, phi, pol, 0)
        heff = reportPtr.GaintoHeight(gain, f, 1.79)
        filter_gain = detectorPtr.GetElectGain_1D_OutZero(f/1e6)
        filter_phase = detectorPtr.GetElectPhase_1D(f/1e6)
        if(np.isnan(heff)):
            heff = 0
        phase = detectorPtr.GetAntPhase_1D(f/1e6, theta, phi, pol)
        gains.append(gain)
        heffs.append(heff*complex(np.cos(np.radians(phase)),np.sin(np.radians(phase))))
        filter_gains.append(filter_gain*complex(np.cos(np.radians(filter_phase)),np.sin(np.radians(filter_phase))))
        # phases.append(phase)
    return np.array(freq),np.array(heffs),np.array(filter_gains)


def deConvolve_antennaAraSim(time, voltage, theta, phi, pol_ant):
    """
    Apply inverse of ARA antenna response
    ----------
    time : array_like
        1D array of times (ns)
    voltage : array_like
        1D array of amplitudes (mV).

    theta, phi, pol_ant: floats
    theta_antenna (radians), phi_antenna (radians), pol_antenna [0:vpol, 1:hpol]
    Returns
    -------
    time : array_like
        1D array of times (ns)
    deDis_wf : array_like
        1D array of amplitudes (mV) of de-convolved waveform.
    """
    import scipy.signal as signal

    sampRate = len(time)/(max(time)-min(time))
    b,a = signal.bessel(4, [0.15,0.4], 'bandpass', analog=False, fs=sampRate)
    fft_v, fft_f, dT = doFFT(time,voltage)
    ff, heffs, filter_gains = getResponseAraSim(theta,phi,fft_f,pol_ant)
    response = heffs
    deDis_wf = np.divide(fft_v,abs(response))
    response = np.divide(response,abs(response))
    deDis_wf = np.divide(deDis_wf,response)
    deDis_wf = np.nan_to_num(deDis_wf)
    revert = doInvFFT(deDis_wf)
    # deDis_wf = signal.lfilter(b, a, revert)
    deDis_wf = revert
    return time, deDis_wf

def weird_division(n, d):
    return n / d if d else 0

def deConvolve(time, voltage, theta, phi, pol_ant):
    """
    Apply inverse of ARA antenna response
    ----------
    time : array_like
        1D array of times (ns)
    voltage : array_like
        1D array of amplitudes (mV).

    theta, phi, pol_ant: floats
    theta_antenna (radians), phi_antenna (radians), pol_antenna [0:vpol, 1:hpol]
    Returns
    -------
    time : array_like
        1D array of times (ns)
    deDis_wf : array_like
        1D array of amplitudes (mV) of de-convolved waveform.
    """
    import scipy.signal as signal
    polarization=np.array([-np.sin(phi),np.cos(phi),-1/np.sin(theta)])
    if(pol_ant == 0):
        ant = ara.VpolAntenna(name="Dummy Vpol", position=(0, 0, 0), power_threshold=0)
        # ant.set_orientation(z_axis=(0, 0, 1), x_axis=(1, 0, 0))#Adding to convert from global coordinates to local antenna coords.
    elif(pol_ant == 1):
        ant = ara.HpolAntenna(name="Dummy Hpol", position=(0, 0, 0), power_threshold=0)
        # ant.set_orientation(z_axis=(0, 0, 1), x_axis=(1, 0, 0))

    sampRate = len(time)/(max(time)-min(time))
    b,a = signal.bessel(4, [0.15,0.4], 'bandpass', analog=False, fs=sampRate)
    fft_v, fft_f, dT = doFFT(time,voltage)
    response_filter = np.array(ant.interpolate_filter(fft_f*1E6))
    dir_res = ant.antenna.directional_response(theta=theta, phi=phi, polarization=polarization)(fft_f*1E6)
    heff = ant.antenna.frequency_response(fft_f*1E6)
    response_antenna = dir_res*heff
    response = response_filter
    # deDis_wf = np.divide(fft_v,abs(response))
    response = np.divide(response,abs(response))
    deDis_wf = np.divide(fft_v,response)
    deDis_wf = np.nan_to_num(deDis_wf)
    revert = doInvFFT(deDis_wf)
    deDis_wf = signal.lfilter(b, a, revert)
    # deDis_wf = revert
    return time, deDis_wf

def PolVectorRecoPower(powerV, powerH, theta, phi):

    R = np.sqrt(powerH/powerV)
    denom = np.sqrt(1+R**2)
    Px = -(np.cos(theta)*np.cos(phi)-R*np.sin(phi))/denom
    Py = -(R*np.cos(phi)+np.cos(theta)*np.sin(phi))/denom
    Pz = np.sin(theta)/denom
    np.set_printoptions(suppress=True)
    # if Peak_V>0:
    #     Px = -Px
    #     Py = -Py
    #     Pz = -Pz
    return np.array([Px,Py,Pz])

def PolVectorRecoPower_signR(powerV, powerH, theta, phi, signR):

    R = np.sqrt(powerH/powerV)*signR
    denom = np.sqrt(1+R**2)
    Px = -(np.cos(theta)*np.cos(phi)-R*np.sin(phi))/denom
    Py = -(R*np.cos(phi)+np.cos(theta)*np.sin(phi))/denom
    Pz = np.sin(theta)/denom
    np.set_printoptions(suppress=True)
    return np.array([Px,Py,Pz])

#Created function for when we already have R calculated.
def PolVectorRecoPower_signR2(R, theta, phi):
    signR = np.sign(R)
    R = R*signR
    denom = np.sqrt(1+R**2)
    Px = -(np.cos(theta)*np.cos(phi)-R*np.sin(phi))/denom
    Py = -(R*np.cos(phi)+np.cos(theta)*np.sin(phi))/denom
    Pz = np.sin(theta)/denom
    np.set_printoptions(suppress=True)
    return np.array([Px,Py,Pz])

def PolVectorReco_debug(Peak_V, Peak_H, theta, phi, sign):
    R = abs(Peak_H/Peak_V)*sign
    denom = np.sqrt(1+R**2)
    Px = -(np.cos(theta)*np.cos(phi)-R*np.sin(phi))/denom
    Py = -(R*np.cos(phi)+np.cos(theta)*np.sin(phi))/denom
    Pz = np.sin(theta)/denom
    np.set_printoptions(suppress=True)
    return np.array([Px,Py,Pz])

def findHighestPeakBin(values):
    if abs(max(values))>=abs(min(values)):
        peakBin = np.argmax(values)
    else:
        peakBin = np.argmin(values)
    return peakBin

def integratePowerWindow(times, values):
    times = np.array(times)
    values = np.array(values)
    dT = times[1]-times[0]
    leftNumBins = int(20/dT)#Number of bins in 20 ns
    rightNumBins = int(60/dT)#Number of bins in 60 ns

    peakBin = findHighestPeakBin(values)#Find bin where peak happens
    lowerEdgeBin = peakBin-leftNumBins
    upperEdgeBin = peakBin+rightNumBins
    if((lowerEdgeBin<0) or (upperEdgeBin<0)):
        return -1
    cutWform = values[lowerEdgeBin:upperEdgeBin]
    cutTimes = times[lowerEdgeBin:upperEdgeBin]
    power = np.sum(cutWform**2)*dT
    return power

def integratePowerWindowSpiceFirstPeak(times, values, maxTime=400):
    times = np.array(times)
    values = np.array(values)
    dT = times[1]-times[0]
    leftNumBins = int(20/dT)#Number of bins in 20 ns
    rightNumBins = int(60/dT)#Number of bins in 60 ns
    peakBin = findHighestPeakBin(values[times<maxTime])#Find bin where peak happens
    lowerEdgeBin = peakBin-leftNumBins
    upperEdgeBin = peakBin+rightNumBins
    if((lowerEdgeBin<0) or (upperEdgeBin<0)):
        return -1
    cutWform = values[lowerEdgeBin:upperEdgeBin]
    cutTimes = times[lowerEdgeBin:upperEdgeBin]
    power = np.sum(cutWform**2)*dT
    return power

def integratePowerWindowSpiceHPol(times, values, startTime, endTime,spiceShift=14.1):
    passTimeCut = True
    times = np.array(times)
    values = np.array(values)
    dT = times[1]-times[0]
    #startTime -= spiceShift
    #endTime -= spiceShift
    #if (startTime < 0):
    #    startTime = 0
    #    endTime = 80
    #lowerEdgeBin = int(startTime/dT)
    #upperEdgeBin = int(endTime/dT)
    #Find number of bins to shift by
    binShift = int(spiceShift/dT)
    lowerEdgeBin = int(startTime - binShift)
    upperEdgeBin = int(endTime - binShift)
    if((lowerEdgeBin<0) or (upperEdgeBin<0)):
        lowerEdgeBin = 0
        upperEdgeBin = int(80/dT)
        passTimeCut = False
        #return -1
    cutWform = values[lowerEdgeBin:upperEdgeBin]
    cutTimes = times[lowerEdgeBin:upperEdgeBin]
    power = np.sum(cutWform**2)*dT
    return power, lowerEdgeBin, upperEdgeBin, passTimeCut

def integratePowerWindowWithTimeCut(times, values, maxTime=400):
    passTimeCut = True
    times = np.array(times)
    values = np.array(values)
    dT = times[1]-times[0]
    leftNumBins = int(20/dT)#Number of bins in 20 ns
    rightNumBins = int(60/dT)#Number of bins in 60 ns
    maxBin = int((maxTime-times[0])/dT) #Find max bin to serve as cutoff for integration of first waveform
    peakBin = findHighestPeakBin(values[:maxBin])#Find bin where peak happens
    lowerEdgeBin = peakBin-leftNumBins
    upperEdgeBin = peakBin+rightNumBins
    if((lowerEdgeBin<0) or (upperEdgeBin<0)):
        lowerEdgeBin = 0
        upperEdgeBin = int(80/dT)
        passTimeCut = False
        #return -1
    if (upperEdgeBin >= len(values)):
        upperEdgeBin = len(values)-1
        lowerEdgeBin = upperEdgeBin - int(80/dT)
        passTimeCut = False
    cutWform = values[lowerEdgeBin:upperEdgeBin]
    cutTimes = times[lowerEdgeBin:upperEdgeBin]
    power = np.sum(cutWform**2)*dT
    return power, lowerEdgeBin, upperEdgeBin, passTimeCut

def integratePowerNoise(times, values):
    times = np.array(times)
    values = np.array(values)
    dT = times[1]-times[0]
    #need to integrate first 80 ns of the waveform
    numBins = int(80/dT) #Trying some debugging - JCF 3/26/2022
    #Integrating over the first 80ns grabs the SpiceCore peak, so we'll integrate over the first 20ns and multiply by 4.
    #numBins = int(20/dT)
    cutWform = values[0:numBins]
    power = np.sum(cutWform**2)*dT
    #return power #Trying some debugging - JCF 3/26/2022
    return power

def integratePowerNoiseFromFullWaveform(times, values):
    times = np.array(times)
    values = np.array(values)
    dT = times[1]-times[0]
    #need to integrate first 80 ns of the waveform
    numBins = int(80/dT) #Trying some debugging - JCF 3/26/2022
    #Integrating over the first 80ns grabs the SpiceCore peak, so we'll integrate over the first 20ns and multiply by 4.
    #numBins = int(20/dT)
    cutWform = values
    power = np.sum(cutWform**2)*dT*numBins/len(cutWform)
    #return power #Trying some debugging - JCF 3/26/2022
    return power

def integratePowerNoiseBasedOnIntegrationWindow(times, values, powerIntegrationStartBin, powerIntegrationEndBin):
    times = np.array(times)
    values = np.array(values)
    dT = times[1]-times[0]
    #Acquire a noise sample over 20ns instead of 80ns.  We'll multiply this power noise by four to match the time window.
    numBins = int(20/dT)
    #If the direct pulse power integration occurs in the first 20ns, we'll get the noise from the last 20ns of the waveform.
    if (numBins > powerIntegrationStartBin):
        cutWform = values[-numBins:]
    #If the direct pulse integration doesn't occur in the first 20ns, we can grab our noise from there.
    else:
        cutWform = values[0:numBins]
    power = np.sum(cutWform**2)*dT
    rms = np.array(cutWform).std()
    #return power noise over 80ns since that matches our power integration window
    return 4*power, rms

def integratePowerNoiseDebug(times, values):
    times = np.array(times)
    values = np.array(values)
    dT = times[1]-times[0]
    ##need to integrate first 80 ns of the waveform
    #numBins = int(80/dT) #Trying some debugging - JCF 3/26/2022
    #Integrating over the first 80ns grabs the SpiceCore peak, so we'll integrate over the first 20ns and multiply by 4.
    numBins = int(20/dT)
    cutWform = values[0:numBins]
    power = np.sum(cutWform**2)*dT
    #return power #Trying some debugging - JCF 3/26/2022
    return 4*power

def integratePowerNoiseEndOfWaveformDebug(times, values):
    times = np.array(times)
    values = np.array(values)
    dT = times[1]-times[0]
    ##need to integrate first 80 ns of the waveform
    #numBins = int(80/dT) #Trying some debugging - JCF 3/26/2022
    #Integrating over the first 80ns grabs the SpiceCore peak, so we'll integrate over the first 20ns and multiply by 4.
    numBins = int(20/dT)
    cutWform = values[-numBins:]
    power = np.sum(cutWform**2)*dT
    #return power #Trying some debugging - JCF 3/26/2022
    return 4*power

def get_view_angles(reportPtr):
	view_angles = {} # return it as a map (which is generally against my rules...)
	for iS in range(4):
		for iA in range(4):
			this_view_angles = []
			view_angle = reportPtr.stations[0].strings[iS].antennas[iA].view_ang
			for v in view_angle:
				this_view_angles.append(v)
			view_angles[(iS*4)+iA] = this_view_angles
	return view_angles

def guess_triggering_solution(eventPtr, reportPtr):
    """
    ​This function tries to estimate what solution -- direct or reflected/refracted --
    likely caused the trigger. What we do is lookup the veiwing angles
    for both solutions, and figure out which which is closer to the Cherenkov cone.
    It turns out that the likely triggering solution is strongly correlated with arrival direction.
    In that if the refracted/reflected ray triggered, the arrival angle of the signal is <85 degrees in zenith.
    And if the direct ray triggered, the angle is >85 degrees in zenith.
    In a real experiment, we would try to measure this arrival angle.
    But this function circumvents doing the incoming wave direction reco.
    If a determination can't be made for some reason, assume a direct solution.
    """
    count_dir = 0 # counter for direct and reflected/refracted rays
    count_ref = 0

    changle_deg = np.rad2deg(eventPtr.Nu_Interaction[0].changle)
    view_angles = get_view_angles(reportPtr)
    for i in range(len(view_angles)):
    	num_sols = len(view_angles[i])
    	if (num_sols != 0 and num_sols !=2):
    		continue
    	elif num_sols==2:
    		dir_angle = changle_deg - np.rad2deg(view_angles[i][0])
    		ref_angle = changle_deg - np.rad2deg(view_angles[i][1])
    		# print("Ant {}, Dir Angle {:.2f}, Ref Angle {:.2f}".format(i, dir_angle, ref_angle))
    		if abs(dir_angle) < abs(ref_angle):
    			count_dir += 1
    		elif abs(ref_angle) < abs(dir_angle):
    			count_ref += 1

    likely_sol = 0 # default to direct
    if count_ref > count_dir:
    	likely_sol = 1
    # print("likely sol {}".format(likely_sol))
    return likely_sol


####Functions created by Justin######
def powerFromWaveform(rawEvent, usefulEvent, vertexReco, ROOT, powerNoiseConstant, gainBalance = False, monteCarlo=False):
    #Import angles
    thetaReco = np.degrees(np.array(vertexReco.reco_arrivalThetas_out))
    phiReco = np.degrees(np.array(vertexReco.reco_arrivalPhis_out))
    
    #Initialize arrays
    powerV = np.zeros(8)
    powerH = np.zeros(8)
    power = np.zeros(16)
    passCutTime = np.zeros(16)
    snrsOut = np.zeros(16)
    passTimeCut = np.ones((16)).astype(bool)
    
    #Initiliaze arrays for power pre noise subtraction and noise.
    powerVPreNoiseSubtraction = np.zeros(8)
    powerVNoiseFromWaveform = np.zeros(8)
    powerHPreNoiseSubtraction = np.zeros(8)
    powerHNoiseFromWaveform = np.zeros(8)
    powerPreNoiseSubtraction = np.zeros(16)
    powerNoiseFromWaveform = np.zeros(16)
    startTime = np.zeros(16)
    endTime = np.zeros(16)
    
    #Check for non-physical angles
    thetaReco[thetaReco<-90] %= 180
    phiReco[phiReco<0]%=360
    
    #Import cutofftime
    try:
        cutoffTime = vertexReco.cutoffTime
    except AttributeError:
        cutoffTime = np.empty(16)
        for ch in range(0,16):
            voltage, time = extractChannelWaveform(usefulEvent, ch)
            cutoffTime[ch] = time[-1]
        
    #Loop over channels
    for ch in range(0,16):
        gr = usefulEvent.getGraphFromRFChan(ch)
        gr = ROOT.FFTtools.getInterpolatedGraph(gr,0.5)
        lenGraph = gr.GetN()    
        time = np.zeros(lenGraph)
        voltage = np.zeros(lenGraph)
        for k in range(0,gr.GetN()):
            time[k] = gr.GetX()[k]
            voltage[k] = gr.GetY()[k]
        if (ch<8):  #VPol case
            #Perform deconvolution
            deConv_t,deConv_v = deConvolve_antenna(time, voltage, np.radians(thetaReco[ch]), np.radians(phiReco[ch]), 0)
            #Perform power integration
            powerPreNoiseSubtraction[ch], startBin, endBin, passTimeCut[ch] = integratePowerWindowWithTimeCut(deConv_t,deConv_v, maxTime=cutoffTime[ch])
            
        else:  #HPol case
            #Perform deconvolution
            deConv_t,deConv_v = deConvolve_antenna(time, voltage, np.radians(thetaReco[ch]), np.radians(phiReco[ch]), 1)
            #Perform power integration
            powerPreNoiseSubtraction[ch], startBin, endBin, passTimeCut[ch] = integratePowerWindowSpiceHPol(deConv_t,deConv_v,startTime[ch-8],endTime[ch-8])
        
        #Calculate noise from waveform
        powerNoiseFromWaveform[ch], rms = integratePowerNoiseBasedOnIntegrationWindow(deConv_t,deConv_v, startBin, endBin)
        # powerNoiseFromWaveform[ch] = integratePowerNoise(deConv_t,deConv_v)
            
        #Find peak voltage
        Peak = findAveragePeak(np.array(deConv_v))
        #Calculate Snr
        snrsOut[ch] = abs(Peak/rms)
        #Perform noise subtraction using channel-based noise from imported distributions.
        power[ch] = powerPreNoiseSubtraction[ch] - powerNoiseConstant[ch]
        
        #Store waveform start and end time for use in HPol integration.
        startTime[ch] = startBin
        endTime[ch] = endBin
        
    if (gainBalance):
        gainCorrection = powerNoiseConstant[:8]/powerNoiseConstant[8:]
        powerPreNoiseSubtraction[8:] *= gainCorrection
        powerNoiseFromWaveform[8:] *= gainCorrection
        power[8:] *= gainCorrection
    #Return the power with and without noise subtraction, the noise from the waveform, and the snrs.
    return power, powerPreNoiseSubtraction, powerNoiseFromWaveform, snrsOut

def noiseFromWaveform(rawEvent, usefulEvent, vertexReco, ROOT, powerNoiseConstant=None, gainBalance = False, monteCarlo=False):
    #Import angles
    thetaReco = np.degrees(np.array(vertexReco.reco_arrivalThetas_out))
    phiReco = np.degrees(np.array(vertexReco.reco_arrivalPhis_out))
    
    #Initialize arrays
    powerV = np.zeros(8)
    powerH = np.zeros(8)
    power = np.zeros(16)
    passCutTime = np.zeros(16)
    snrsOut = np.zeros(16)
    passTimeCut = np.ones((16)).astype(bool)
    
    #Initiliaze arrays for power pre noise subtraction and noise.
    powerVPreNoiseSubtraction = np.zeros(8)
    powerVNoiseFromWaveform = np.zeros(8)
    powerHPreNoiseSubtraction = np.zeros(8)
    powerHNoiseFromWaveform = np.zeros(8)
    powerPreNoiseSubtraction = np.zeros(16)
    powerNoiseFromWaveform = np.zeros(16)
    startTime = np.zeros(16)
    endTime = np.zeros(16)
    
    #Check for non-physical angles
    thetaReco[thetaReco<-90] %= 180
    phiReco[phiReco<0]%=360
    
    #Import cutofftime
    try:
        cutoffTime = vertexReco.cutoffTime
    except AttributeError:
        cutoffTime = np.empty(16)
        for ch in range(0,16):
            voltage, time = extractChannelWaveform(usefulEvent, ch)
            cutoffTime[ch] = time[-1]
        
    #Loop over channels
    for ch in range(0,16):
        gr = usefulEvent.getGraphFromRFChan(ch)
        gr = ROOT.FFTtools.getInterpolatedGraph(gr,0.5)
        lenGraph = gr.GetN()    
        time = np.zeros(lenGraph)
        voltage = np.zeros(lenGraph)
        for k in range(0,gr.GetN()):
            time[k] = gr.GetX()[k]
            voltage[k] = gr.GetY()[k]
        if (ch<8):  #VPol case
            #Perform deconvolution
            deConv_t,deConv_v = deConvolve_antenna(time, voltage, np.radians(thetaReco[ch]), np.radians(phiReco[ch]), 0)
            #Perform power integration
            powerPreNoiseSubtraction[ch], startBin, endBin, passTimeCut[ch] = integratePowerWindowWithTimeCut(deConv_t,deConv_v, maxTime=cutoffTime[ch])
            
        else:  #HPol case
            #Perform deconvolution
            deConv_t,deConv_v = deConvolve_antenna(time, voltage, np.radians(thetaReco[ch]), np.radians(phiReco[ch]), 1)
            #Perform power integration
            powerPreNoiseSubtraction[ch], startBin, endBin, passTimeCut[ch] = integratePowerWindowSpiceHPol(deConv_t,deConv_v,startTime[ch-8],endTime[ch-8])
        
        #Calculate noise from waveform
        powerNoiseFromWaveform[ch], rms = integratePowerNoiseBasedOnIntegrationWindow(deConv_t,deConv_v, startBin, endBin)
        # powerNoiseFromWaveform[ch] = integratePowerNoise(deConv_t,deConv_v)
            
        #Find peak voltage
        Peak = findAveragePeak(np.array(deConv_v))
        #Calculate Snr
        snrsOut[ch] = abs(Peak/rms)
        #Perform noise subtraction using channel-based noise from imported distributions.
        if (powerNoiseConstant is None):
            power[ch] = powerPreNoiseSubtraction[ch]
        else:
            power[ch] = powerPreNoiseSubtraction[ch] - powerNoiseConstant[ch]
            
            if (gainBalance):
                gainCorrection = powerNoiseConstant[:8]/powerNoiseConstant[8:]
                powerPreNoiseSubtraction[8:] *= gainCorrection
                powerNoiseFromWaveform[8:] *= gainCorrection
                power[8:] *= gainCorrection
        
        #Store waveform start and end time for use in HPol integration.
        startTime[ch] = startBin
        endTime[ch] = endBin
        

    #Return the power with and without noise subtraction, the noise from the waveform, and the snrs.
    return power, powerPreNoiseSubtraction, powerNoiseFromWaveform, snrsOut
        
def calculatePsiAndR(powerArray):
    powerH = powerArray[8:]
    powerV = powerArray[:8]
    powerSign = np.ones(8)
    for pair in range(8):
        if (powerH[pair] < 0 or powerV[pair] < 0):
            powerSign[pair] = -1
    R = np.sqrt(abs(powerH/powerV))*powerSign
    Psi = np.degrees(np.arctan(R))
    return R, Psi

def calculateMeanNoise(eventList):
    allEnvelopes = np.zeros(16)
    for evt in eventList:
        #Import RF Event and subtract noise.
        eventTree.GetEntry(evt)
        usefulEvent = ROOT.UsefulAtriStationEvent(rawEvent,ROOT.AraCalType.kLatestCalib)
        softEventNumber = usefulEvent.eventNumber
        noiseMean = np.empty(16)
        for ch in range(0,16):
            gr = usefulEvent.getGraphFromRFChan(ch)
            gr = ROOT.FFTtools.getInterpolatedGraph(gr,0.5)
            lenGraph = gr.GetN()
            t = np.zeros(lenGraph)
            v = np.zeros(lenGraph)
            for k in range(0,gr.GetN()):
                t[k] = gr.GetX()[k]
                v[k] = gr.GetY()[k]
            analytical_signal_soft = hilbert(v)  #Real component is the initial data.  Imaginary component is the Hilbert transform
            amplitude_envelope_soft = np.abs(analytical_signal_soft) #Magnitude is the envelope of the function
            noiseMean[ch] = amplitude_envelope_soft.mean()
            allEnvelopes[ch] += amplitude_envelope_soft.mean()    
    meanEnvelope = allEnvelopes/len(eventList)
    
    return meanEnvelope

def powerFromWaveformSubtractHilbertNoise(rawEvent, usefulEvent, vertexReco, ROOT, noiseEnvelope=None, noiseRms=None, gainBalance=False, gainCorrection=None):
    from scipy.signal import hilbert
    #TODO:  Need to write gainBalance for finding noise from the waveform.
    
    #Import angles
    thetaReco = np.degrees(np.array(vertexReco.reco_arrivalThetas_out))
    phiReco = np.degrees(np.array(vertexReco.reco_arrivalPhis_out))
    
    #Initialize arrays
    powerV = np.zeros(8)
    powerH = np.zeros(8)
    power = np.zeros(16)
    passCutTime = np.zeros(16)
    snrsOut = np.zeros(16)
    passTimeCut = np.ones((16)).astype(bool)   
    
    #Initiliaze arrays for power pre noise subtraction and noise.
    powerVPreNoiseSubtraction = np.zeros(8)
    powerVNoiseFromWaveform = np.zeros(8)
    powerHPreNoiseSubtraction = np.zeros(8)
    powerHNoiseFromWaveform = np.zeros(8)
    powerPreNoiseSubtraction = np.zeros(16)
    powerNoiseFromWaveform = np.zeros(16)
    startTime = np.zeros(16)
    endTime = np.zeros(16)
    
    #Check for non-physical angles
    thetaReco[thetaReco<-90] %= 180
    phiReco[phiReco<0]%=360
    
    #Import cutofftime
    try:
        cutoffTime = vertexReco.cutoffTime
    except AttributeError:
        cutoffTime = np.empty(16)
        for ch in range(0,16):
            voltage, time = extractChannelWaveform(usefulEvent, ch)
            cutoffTime[ch] = time[-1]
            
        
    #Loop over channels
    for ch in range(0,16):
        gr = usefulEvent.getGraphFromRFChan(ch)
        gr = ROOT.FFTtools.getInterpolatedGraph(gr,0.5)
        lenGraph = gr.GetN()    
        time = np.zeros(lenGraph)
        voltage = np.zeros(lenGraph)
        for k in range(0,gr.GetN()):
            time[k] = gr.GetX()[k]
            voltage[k] = gr.GetY()[k]
            
        #Calculate SNR before doing the Hilbert enveloping.
        if (ch<8):  #VPol case
            #Perform deconvolution
            deConv_t,deConv_v = deConvolve_antenna(time, voltage, np.radians(thetaReco[ch]), np.radians(phiReco[ch]), 0)
            #Perform power integration
            power[ch], startBin, endBin, passTimeCut[ch] = integratePowerWindowWithTimeCut(deConv_t,deConv_v, maxTime=cutoffTime[ch])
            
        else:  #HPol case
            #Perform deconvolution
            deConv_t,deConv_v = deConvolve_antenna(time, voltage, np.radians(thetaReco[ch]), np.radians(phiReco[ch]), 1)
            #Perform power integration
            power[ch], startBin, endBin, passTimeCut[ch] = integratePowerWindowSpiceHPol(deConv_t,deConv_v,startTime[ch-8],endTime[ch-8])
    
        #Find peak voltage
        Peak = findAveragePeak(np.array(deConv_v))
            
        #Calculate Hilbert envelope of waveform
        analytical_signal_rf = hilbert(voltage)  #Real component is the initial data.  Imaginary component is the Hilbert transform
        amplitude_envelope_rf = np.abs(analytical_signal_rf) #Magnitude is the envelope of the function
        
        #Perform noise subtraction of Hilbert Envelope
        if (noiseEnvelope is None):
            envelopeAfterSubtraction = amplitude_envelope_rf
        else:
            envelopeAfterSubtraction = amplitude_envelope_rf - noiseEnvelope[ch]
        envelopePeak = max(envelopeAfterSubtraction)
        #Find noise reduction ratio
        noiseReductionRatio = envelopeAfterSubtraction/amplitude_envelope_rf
        #Apply noise reduction ratio to waveform
        voltage *= noiseReductionRatio        
        
        #Perform deconvolution after noise subtraction.
        if (ch<8):  #VPol case
            #Perform deconvolution
            deConv_t,deConv_v = deConvolve_antenna(time, voltage, np.radians(thetaReco[ch]), np.radians(phiReco[ch]), 0)
            #Perform power integration
            power[ch], startBin, endBin, passTimeCut[ch] = integratePowerWindowWithTimeCut(deConv_t,deConv_v, maxTime=cutoffTime[ch])
            
        else:  #HPol case
            #Perform deconvolution
            deConv_t,deConv_v = deConvolve_antenna(time, voltage, np.radians(thetaReco[ch]), np.radians(phiReco[ch]), 1)
            #Perform power integration
            power[ch], startBin, endBin, passTimeCut[ch] = integratePowerWindowSpiceHPol(deConv_t,deConv_v,startTime[ch-8],endTime[ch-8])
        #Calculate Snr
        snrsOut[ch] = envelopePeak/noiseRms[ch]
    
        #Store waveform start and end time for use in HPol integration.
        startTime[ch] = startBin
        endTime[ch] = endBin
        
    #Gain balance needs to be discussed with Amy before we implement it.
    if (gainBalance):
        gainCorrectionFactor = gainCorrection[:8]/gainCorrection[8:]
        power[8:] *= gainCorrectionFactor
    
    #return power and snr
    return power, snrsOut


def findMeanNoise(eventList, eventTree, rawEvent, ROOT, waveformSampleNs = None):
    from scipy.signal import hilbert
    allEnvelopes = np.zeros(16)
    numEntries = np.zeros(16)
    rmsSquared = np.zeros(16)
    for evt in eventList:
        #Import RF Event and subtract noise.
        eventTree.GetEntry(evt)
        usefulEvent = ROOT.UsefulAtriStationEvent(rawEvent,ROOT.AraCalType.kLatestCalib)
        softEventNumber = usefulEvent.eventNumber
        noiseMean = np.empty(16)
        
        for ch in range(0,16):
            gr = usefulEvent.getGraphFromRFChan(ch)
            gr = ROOT.FFTtools.getInterpolatedGraph(gr,0.5)
            lenGraph = gr.GetN()
            numEntries[ch] += lenGraph
            #t = []
            #v = []
            t = np.zeros(lenGraph)
            v = np.zeros(lenGraph)
            for k in range(0,gr.GetN()):
                #t.append(gr.GetX()[k])
                #v.append(gr.GetY()[k])
                t[k] = gr.GetX()[k]
                v[k] = gr.GetY()[k]
            if (waveformSampleNs is not None):
                t = t[:2*waveformSampleNs]
                v = v[:2*waveformSampleNs]
            analytical_signal_soft = hilbert(v)  #Real component is the initial data.  Imaginary component is the Hilbert transform
            amplitude_envelope_soft = np.abs(analytical_signal_soft) #Magnitude is the envelope of the function
            noiseMean[ch] = amplitude_envelope_soft.mean()
            rmsSquared[ch] += (v.std())**2
            allEnvelopes[ch] += amplitude_envelope_soft.mean()
    
    meanEnvelope = allEnvelopes/len(eventList)
    # meanEnvelope = allEnvelopes/numEntries
    rms = np.sqrt(rmsSquared/len(eventList))
    
    if (len(eventList)==0):
        meanEnvelope[:] = 0
        rms[:] = 0
    
    return meanEnvelope, rms

#adding noise calculation method for simulated pulser events by simply grabbing noise from beginning of waveform.
def findMeanNoiseFromWaveform(eventList, eventTree, usefulEvent, ROOT, sampleTime=40):
    from scipy.signal import hilbert
    allEnvelopes = np.zeros(16)
    numEntries = np.zeros(16)
    rmsSquared = np.zeros(16)
    
    for evt in eventList:
        #Import RF Event and subtract noise.
        eventTree.GetEntry(evt)
        # usefulEvent = ROOT.UsefulAtriStationEvent(rawEvent,ROOT.AraCalType.kLatestCalib)
        # softEventNumber = usefulEvent.eventNumber
        noiseMean = np.empty(16)
        
        for ch in range(0,16):
            gr = usefulEvent.getGraphFromRFChan(ch)
            gr = ROOT.FFTtools.getInterpolatedGraph(gr,0.5)
            # lenGraph = gr.GetN()
            lenGraph = sampleTime*2 #Grabs a sampleTime's (default 40ns) worth of data from beginning of waveform.
            numEntries[ch] += lenGraph
            #t = []
            #v = []
            t = np.zeros(lenGraph)
            v = np.zeros(lenGraph)
            for k in range(0,lenGraph):
                #t.append(gr.GetX()[k])
                #v.append(gr.GetY()[k])
                t[k] = gr.GetX()[k]
                v[k] = gr.GetY()[k]
            analytical_signal_soft = hilbert(v)  #Real component is the initial data.  Imaginary component is the Hilbert transform
            amplitude_envelope_soft = np.abs(analytical_signal_soft) #Magnitude is the envelope of the function
            noiseMean[ch] = amplitude_envelope_soft.mean()
            rmsSquared[ch] += (v.std())**2
            allEnvelopes[ch] += amplitude_envelope_soft.mean()
    
    meanEnvelope = allEnvelopes/len(eventList)
    # meanEnvelope = allEnvelopes/numEntries
    rms = np.sqrt(rmsSquared/len(eventList))
    
    if (len(eventList)==0):
        meanEnvelope[:] = 0
        rms[:] = 0
    
    return meanEnvelope, rms

def powerFromSoftTriggerNoDeconvolution(rawEvent, usefulEvent, ROOT):
    
    #Initialize arrays
    power = np.zeros(16)
    passCutTime = np.zeros(16)
    snrsOut = np.zeros(16)
    passTimeCut = np.ones((16)).astype(bool)   
    
    #Initiliaze arrays for power pre noise subtraction and noise.
    powerVPreNoiseSubtraction = np.zeros(8)
    powerVNoiseFromWaveform = np.zeros(8)
    powerHPreNoiseSubtraction = np.zeros(8)
    powerHNoiseFromWaveform = np.zeros(8)
    powerPreNoiseSubtraction = np.zeros(16)
    powerNoiseFromWaveform = np.zeros(16)
    startTime = np.zeros(16)
    endTime = np.zeros(16)
        
    #Loop over channels
    for ch in range(0,16):
        gr = usefulEvent.getGraphFromRFChan(ch)
        gr = ROOT.FFTtools.getInterpolatedGraph(gr,0.5)
        lenGraph = gr.GetN()    
        time = np.zeros(lenGraph)
        voltage = np.zeros(lenGraph)
        for k in range(0,gr.GetN()):
            time[k] = gr.GetX()[k]
            voltage[k] = gr.GetY()[k]
            
        #Calculate power.
        power[ch] = integratePowerNoise(time,voltage)
    #return power
    return power

# def peakHilbert(rawEvent, usefulEvent, vertexReco, noiseEnvelope, noiseRms, gainBalance=False, gainCorrection=None, tolerance=10, timeShift=14.1, debug=False, deconvolution=True, solution='D'):
def peakHilbert(usefulEvent, vertexReco, noiseEnvelope, noiseRms, gainBalance=False, gainCorrection=None, tolerance=10, timeShift=14.1, debug=False, deconvolution=True, solution='D', station=2, configuration=None):
    #TODO: Update this to find peak from HPol when HPol is larger than VPol, or by user command (some kind off flag)
    import ROOT
    from scipy.signal import hilbert
    #TODO:  Need to write gainBalance for finding noise from the waveform.
    
    #Import angles
    thetaReco = np.degrees(np.array(vertexReco.reco_arrivalThetas_out))
    phiReco = np.degrees(np.array(vertexReco.reco_arrivalPhis_out))
    
    #Initialize arrays
    powerV = np.zeros(8)
    powerH = np.zeros(8)
    power = np.zeros(16)
    passCutTime = np.zeros(16)
    snrsOut = np.zeros(16)
    passTimeCut = np.ones((16)).astype(bool)   
    
    #Initiliaze arrays for power pre noise subtraction and noise.
    powerVPreNoiseSubtraction = np.zeros(8)
    powerVNoiseFromWaveform = np.zeros(8)
    powerHPreNoiseSubtraction = np.zeros(8)
    powerHNoiseFromWaveform = np.zeros(8)
    powerPreNoiseSubtraction = np.zeros(16)
    powerNoiseFromWaveform = np.zeros(16)
    startTime = np.zeros(16)
    endTime = np.zeros(16)
    peakTime = np.zeros(16)
    hilbertPeak = np.zeros(16)
    peakWindows = np.zeros((8,2))
    
    #Check for non-physical angles
    thetaReco[thetaReco<-90] %= 180
    phiReco[phiReco<0]%=360
    
    #Import cutofftime
    try:
        cutoffTime = vertexReco.cutoffTime
    except AttributeError:
        cutoffTime = np.empty(16)
        for ch in range(0,16):
            voltage, time = extractChannelWaveform(usefulEvent, ch)
            cutoffTime[ch] = time[-1]
    #Debugging and apply gain balance before deconvolution - JCF 6/16/2023
    if (gainBalance):
        gainCorrectionFactor = gainCorrection[:8]/gainCorrection[8:]
    #End debugging    
    
    #TODO: Have this loop over channel pairs instead, so that we may find if VPol or HPol is stronger.    
    for channelPair in range(8):
        voltageVpol, timeVpol = extractChannelWaveform(usefulEvent, channelPair)
        voltageHpol, timeHpol = extractChannelWaveform(usefulEvent, channelPair+8)
        print("voltageVpol = " +str(voltageVpol))
        print("voltageHpol = " +str(voltageHpol))
        
        #Debugging and apply gain balance before deconvolution - JCF 6/16/2023
        if (gainBalance):
            voltageHpol*=gainCorrectionFactor[channelPair]
        #End debugging    
        
        #Calculate Hilbert envelope of waveforms
        amplitude_envelope_rf_Vpol = getHilbertEnvelopeSingleChannel(voltageVpol)
        amplitude_envelope_rf_Hpol = getHilbertEnvelopeSingleChannel(voltageHpol)
        
        #Perform noise subtraction of Hilbert Envelopes
        if (noiseEnvelope is not None):
            envelopeAfterSubtractionVpol = amplitude_envelope_rf_Vpol - noiseEnvelope[channelPair]
            envelopeAfterSubtractionHpol = amplitude_envelope_rf_Hpol - noiseEnvelope[channelPair+8]
        else:
            envelopeAfterSubtractionVpol = amplitude_envelope_rf_Vpol
            envelopeAfterSubtractionHpol = amplitude_envelope_rf_Hpol           
        
        #Find noise reduction ratio
        noiseReductionRatioVpol = envelopeAfterSubtractionVpol/amplitude_envelope_rf_Vpol
        noiseReductionRatioHpol = envelopeAfterSubtractionHpol/amplitude_envelope_rf_Hpol
        
        #Apply noise reduction ratio to waveform
        voltageVpol *= noiseReductionRatioVpol
        voltageHpol *= noiseReductionRatioHpol
        
        #Deconvolve waveforms, if desired
        if (deconvolution):
            # Debugging - JCF 6/8/2023
            deConv_t_Vpol,deConv_v_Vpol = deConvolve_antenna(timeVpol, voltageVpol, np.radians(thetaReco[channelPair]), np.radians(phiReco[channelPair]), 0, station=station, channel=channelPair, configuration=configuration)
            deConv_t_Hpol,deConv_v_Hpol = deConvolve_antenna(timeHpol, voltageHpol, np.radians(thetaReco[channelPair+8]), np.radians(phiReco[channelPair+8]), 1, station=station, channel=channelPair+8, configuration=configuration)
            # deConv_t_Vpol,deConv_v_Vpol = deConvolve_antennaAraSim(timeVpol, voltageVpol, np.radians(thetaReco[channelPair]), np.radians(phiReco[channelPair]), 0)
            # deConv_t_Hpol,deConv_v_Hpol = deConvolve_antennaAraSim(timeHpol, voltageHpol, np.radians(thetaReco[channelPair+8]), np.radians(phiReco[channelPair+8]), 1)
            # End debugging
        else:
            deConv_t_Vpol, deConv_v_Vpol = timeVpol, voltageVpol
            deConv_t_Hpol, deConv_v_Hpol = timeHpol, voltageHpol
        
        #Find which channel is stronger.
        peakVpol = max(envelopeAfterSubtractionVpol)
        peakHpol = max(envelopeAfterSubtractionHpol)
        
        print(peakVpol)
        print(peakHpol)
        
        print(deConv_v_Vpol)
        print(deConv_v_Hpol)        
        
        if (peakVpol > peakHpol):
            primaryVoltage = deConv_v_Vpol
            secondaryVoltage = deConv_v_Hpol
            primaryTime = deConv_t_Vpol
            secondaryTime = deConv_t_Hpol
            cutoffTimePair = np.array([cutoffTime[channelPair],cutoffTime[channelPair+8]])
            primaryPolarization = 0
            hilbertPeak[channelPair], hilbertPeak[channelPair+8], peakTime[channelPair], peakTime[channelPair+8] = findHilbertForLargerChannel(primaryVoltage, secondaryVoltage, primaryTime, secondaryTime, cutoffTimePair, primaryPolarization, solution=solution, tolerance=tolerance, timeShift=timeShift)
        else:
            primaryVoltage = deConv_v_Hpol
            secondaryVoltage = deConv_v_Vpol
            primaryTime = deConv_t_Hpol
            secondaryTime = deConv_t_Vpol
            cutoffTimePair = np.array([cutoffTime[channelPair+8],cutoffTime[channelPair]])
            primaryPolarization = 1
            hilbertPeak[channelPair+8], hilbertPeak[channelPair], peakTime[channelPair+8], peakTime[channelPair] = findHilbertForLargerChannel(primaryVoltage, secondaryVoltage, primaryTime, secondaryTime, cutoffTimePair, primaryPolarization, solution=solution, tolerance=tolerance, timeShift=timeShift)
            
    #Calculate Snr.
    snrsOut = hilbertPeak/noiseRms        
        
    #Gain balance needs to be discussed with Amy before we implement it.
    # #Debugging and apply gain balance before deconvolution - JCF 6/16/2023
    # if (gainBalance):
    #     gainCorrectionFactor = gainCorrection[:8]/gainCorrection[8:]
    #     hilbertPeak[8:] *= gainCorrectionFactor
    # #End debugging
    
    #return power and snr
    # if (debug):
    #     return hilbertPeak, peakTime, snrsOut, peakWindows
    # else:
    return hilbertPeak, peakTime, snrsOut
    
def findHilbertForLargerChannel(primaryVoltage, secondaryVoltage, primaryTime, secondaryTime, cutoffTimePair, primaryPolarization, solution='single', tolerance=10, timeShift=0):
    #This function is for use in peakHilbert, which finds the hilbert peaks based on whether the Hpol or Vpol channel has a stronger signal.
    #The stronger channel serves as the primary voltage in this function, while the weaker partner channel is the secondary.
    
    #Calculate hilbert peak of primary channel
    amplitude_envelope_rf_Primary = getHilbertEnvelopeSingleChannel(primaryVoltage)
    if solution in ["D", "d", "Direct", "direct", 0]:
        amplitude_envelope_rf_Primary = amplitude_envelope_rf_Primary[primaryTime<cutoffTimePair[0]]
        primaryTime = primaryTime[primaryTime<cutoffTimePair[0]]
    elif solution in ['R', 'r', 'Reflected', 'reflected', 'Refracted', 'refracted', 1]:
        amplitude_envelope_rf_Primary = amplitude_envelope_rf_Primary[primaryTime>cutoffTimePair[0]]
        primaryTime = primaryTime[primaryTime>cutoffTimePair[0]]
    elif solution in ['single', 's', 'noise', 'calpulser']:
        amplitude_envelope_rf_Primary = amplitude_envelope_rf_Primary
        primaryTime = primaryTime
    else:
        print("Solution choice unknown.  Must choose either direct or reflected/refracted solution.  Exiting.")
        sys.exit()        
    envelopePeakPrimary = max(amplitude_envelope_rf_Primary)
    peakTimePrimary = primaryTime[np.where(amplitude_envelope_rf_Primary == envelopePeakPrimary)]
    
    #Calculate Hilbert peak of secondary channel
    amplitude_envelope_rf_Secondary = getHilbertEnvelopeSingleChannel(secondaryVoltage)
    
    if (np.isscalar(tolerance)):
        minTime = peakTimePrimary-tolerance-timeShift*(1-2*primaryPolarization) 
        maxTime = peakTimePrimary+tolerance-timeShift*(1-2*primaryPolarization)
        #The (1-2*primaryPolarization) term is a trick for changing the sign based on the primary polarization.  If Vpol (primaryPolarization = 0), it subtracts the timeshift.  If Hpol (primaryPolarization = 1), it adds the timeshift.
    else:
        minTime = peakTimePrimary-tolerance[0]-timeShift*(1-2*primaryPolarization)
        maxTime = peakTimePrimary+tolerance[1]-timeShift*(1-2*primaryPolarization)
        
    
    print("primaryVoltage = " + str(primaryVoltage))
    print("amplitude_envelope_rf_Primary = " + str(amplitude_envelope_rf_Primary))
    print("envelopePeakPrimary = " + str(envelopePeakPrimary))
    print("secondaryTime = " + str(secondaryTime))
    print("peakTimePrimary = " + str(peakTimePrimary))    
    print("tolerance = " + str(tolerance))  
    print("timeShift = " + str(timeShift))  
    print("primaryPolarization = " + str(primaryPolarization))  
    print("minTime = " + str(minTime))
    print("maxTime = " + str(maxTime))
    peakWindowSecondary = (secondaryTime>=minTime)*(secondaryTime<=maxTime)
    amplitude_envelope_rf_Secondary = amplitude_envelope_rf_Secondary[peakWindowSecondary]
    
    try:
        envelopePeakSecondary = max(amplitude_envelope_rf_Secondary)
        secondaryTime = secondaryTime[peakWindowSecondary]
        peakTimeSecondary = secondaryTime[np.where(amplitude_envelope_rf_Secondary == envelopePeakSecondary)]
        if not np.isscalar(peakTimeSecondary):
            peakTimeSecondary = np.mean(peakTimeSecondary)
    except ValueError as ve:
        print("Secondary envelope empty.  Setting envelopePeak to zero.")
        envelopePeakSecondary=0
        peakTimeSecondary = 0
        
    # print(min(secondaryTime))
    # print(max(secondaryTime))
    
    # secondaryTime = secondaryTime[peakWindowSecondary]
    # peakTimeSecondary = secondaryTime[np.where(amplitude_envelope_rf_Secondary == envelopePeakSecondary)]
    
    # print(envelopePeakPrimary)
    # print(envelopePeakSecondary)
    # print(peakTimePrimary)
    # print(peakTimeSecondary)
    # print(minTime)
    # print(maxTime)
    
    return envelopePeakPrimary, envelopePeakSecondary, peakTimePrimary, peakTimeSecondary

def extractChannelWaveform(usefulEvent, ch, sampleRate=0.5):
    #Takes ARA useful event, extracts waveform info for the selected channel, and interpolates it to new sample rate (default 0.5 ns).
    import ROOT
    from ROOT import gInterpreter, gSystem
    from ROOT import TChain, TSelector, TTree
    gr = usefulEvent.getGraphFromRFChan(ch)
    gr = ROOT.FFTtools.getInterpolatedGraph(gr,sampleRate)
    lenGraph = gr.GetN()    
    time = np.zeros(lenGraph)
    voltage = np.zeros(lenGraph)
    for k in range(0,gr.GetN()):
        time[k] = gr.GetX()[k]
        voltage[k] = gr.GetY()[k]
    return voltage, time
                
def getHilbertEnvelopeSingleChannel(voltage):
    from scipy.signal import hilbert
    hilbertEnvelope = np.abs(hilbert(voltage)) #Magnitude is the envelope of the function
    return hilbertEnvelope
    
def getHilbertEnvelopesAllChannels(usefulEvent, sampleWindow=400):
    import ROOT
    from scipy.signal import hilbert
    hilbertEnvelopes = np.zeros((16,sampleWindow))
    for ch in range(16):
        voltage, time = extractChannelWaveform(usefulEvent,ch)
        channelEnvelope = getHilbertEnvelopeSingleChannel(voltage)
        hilbertEnvelopes[ch,:] = channelEnvelope[:sampleWindow]
    return hilbertEnvelopes

# def plotRawWaveform(usefulEvent, timeMarkers=None):
#     import ROOT
#     import matplotlib.pyplot as plt
#     import numpy
#     fig, ax = plt.subplots(nrows=4,ncols=4, figsize=(12,12),sharex=True, sharey=True)
#     eventNumber = usefulEvent.eventNumber
#     for ch in range(16):
#         voltage, time = extractChannelWaveform(usefulEvent, ch)
#         ax[int(ch/4), ch%4].set_title("Channel " + str(ch))
#         ax[int(ch/4), ch%4].plot(time,voltage, alpha=0.7)
#     plt.suptitle("Event Number: " +str(eventNumber))
#     if (timeMarkers is not None):
#         if (numpy.isscalar(timeMarkers)):
#             for ch in range(16):
#                 ax[int(ch/4), ch%4].axvline(timeMarkers, color='black', linestyle='--')
#         elif (len(timeMarkers)==16):
#             for ch in range(16):
#                 ax[int(ch/4), ch%4].axvline(timeMarkers[ch], color='black', linestyle='--')

def plotRawWaveform(usefulEvent, timeMarkers=None, runNumber=None, station = None, channelPair=None):
    import ROOT
    import matplotlib.pyplot as plt
    import numpy
    if (channelPair is None):
        fig, ax = plt.subplots(nrows=4,ncols=4, figsize=(12,12),sharex=True, sharey=True)
        channels = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
    elif (np.isscalar(channelPair)):
        fig, ax = plt.subplots(nrows=2,ncols=1, figsize=(12,12),sharex=True, sharey=True)
        channels = [channelPair, channelPair+8]
    fig.supxlabel("Time [ns]")
    fig.supylabel("Voltage [mV]")
    eventNumber = usefulEvent.eventNumber
    for ch, axes in zip(channels, ax.flatten()):
        voltage, time = extractChannelWaveform(usefulEvent, ch)
        axes.set_title("Channel " + str(ch))
        axes.plot(time,voltage, alpha=0.7)    
    # for ch in range(16):
    #     voltage, time = extractChannelWaveform(usefulEvent, ch)
    #     ax[int(ch/4), ch%4].set_title("Channel " + str(ch))
    #     ax[int(ch/4), ch%4].plot(time,voltage, alpha=0.7)
    if (station is not None):
        stationTitle = "Station: " + str(station) + " | "
    else:
        stationTitle = ''
    if (runNumber is not None):
        runTitle = "Run Number: " + str(runNumber) + " | "
    else:
        runTitle = ''
    plt.suptitle(stationTitle + runTitle + "Event Number: " +str(eventNumber))
   
    if (timeMarkers is not None):
        if (numpy.isscalar(timeMarkers)):
            for ch, axes in zip(channels, ax.flatten()):
                axes.axvline(timeMarkers, color='black', linestyle='--')
        elif (len(timeMarkers)==16):
            for ch, axes in zip(channels, ax.flatten()):
                axes.axvline(timeMarkers[ch], color='black', linestyle='--')    
        
        
        # if (numpy.isscalar(timeMarkers)):
        #     for ch in range(16):
        #         ax[int(ch/4), ch%4].axvline(timeMarkers, color='black', linestyle='--')
        # elif (len(timeMarkers)==16):
        #     for ch in range(16):
        #         ax[int(ch/4), ch%4].axvline(timeMarkers[ch], color='black', linestyle='--')
    return fig, ax
                
def plotDeconvolvedWaveform(usefulEvent, vertexReco, timeMarkers=None, runNumber=None, station = None, Hilbert=False, channelPair=None, deconvolution=True, extraInfo = '', showPeaks = False, tolerance=10, rawEvent=None, noiseEnvelope=None, noiseRms=None, noiseReduction=False, configuration=None):
    import ROOT
    import matplotlib.pyplot as plt
    import numpy
    # from scipy.signal import hilbert
    if (channelPair is None):
        fig, ax = plt.subplots(nrows=4,ncols=4, figsize=(12,12),sharex=True, sharey=True)
        channels = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
    elif (np.isscalar(channelPair)):
        fig, ax = plt.subplots(nrows=2,ncols=1, figsize=(12,12),sharex=True, sharey=True)
        channels = [channelPair, channelPair+8]
    fig.supxlabel("Time [ns]")
    fig.supylabel("Voltage [mV]")
    eventNumber = usefulEvent.eventNumber
    
    if (showPeaks and rawEvent is not None):
        peaks, peakTimes, snrs, peakWindows = peakHilbert(rawEvent, usefulEvent, vertexReco, noiseEnvelope, noiseRms, debug=True, tolerance=tolerance, deconvolution=deconvolution)
    
    #Import angles
    thetaReco = np.degrees(np.array(vertexReco.reco_arrivalThetas_out))
    phiReco = np.degrees(np.array(vertexReco.reco_arrivalPhis_out))
    
    for ch, axes in zip(channels, ax.flatten()):
        voltage, time = extractChannelWaveform(usefulEvent, ch)
        if (noiseReduction):    
            amplitude_envelope_rf = getHilbertEnvelopeSingleChannel(voltage)        
            #Perform noise subtraction of Hilbert Envelope
            envelopeAfterSubtraction = amplitude_envelope_rf - noiseEnvelope[ch]
            #Isolate envelope to first peak using the cutoff time
            # envelopeAfterSubtraction = envelopeAfterSubtraction[time<cutoffTime[ch]]
            # envelopePeak = max(envelopeAfterSubtraction)
            #Find noise reduction ratio
            noiseReductionRatio = envelopeAfterSubtraction/amplitude_envelope_rf
            #Apply noise reduction ratio to waveform
            voltage *= noiseReductionRatio
        
        # if (ch<8):  #VPol case
        #     #Perform deconvolution
        #     deConv_t,deConv_v = deConvolve_antenna(time, voltage, np.radians(thetaReco[ch]), np.radians(phiReco[ch]), 0)
        # else:  #HPol case
        #     #Perform deconvolution
        #     deConv_t,deConv_v = deConvolve_antenna(time, voltage, np.radians(thetaReco[ch]), np.radians(phiReco[ch]), 1)
        if (ch<8):  #VPol case
            #Perform deconvolution
            if (deconvolution):
                    #Perform deconvolution
                    deConv_t,deConv_v = deConvolve_antenna(time, voltage, np.radians(thetaReco[ch]), np.radians(phiReco[ch]), 0, station=station, channel=ch, configuration=configuration)
            else:
                deConv_t = time
                deConv_v = voltage
        else:  #HPol case
            #Perform deconvolution
            if (deconvolution):
                    #Perform deconvolution
                    deConv_t,deConv_v = deConvolve_antenna(time, voltage, np.radians(thetaReco[ch]), np.radians(phiReco[ch]), 1, station=station, channel=ch, configuration=configuration)
            else:
                deConv_t = time
                deConv_v = voltage

            
        axes.set_title("Channel " + str(ch))
        axes.plot(deConv_t,deConv_v, alpha=0.7, label='Waveform')
        if (showPeaks and rawEvent is not None):
            axes.axvline(peakTimes[ch], linestyle='--', color='black')
            axes.axhline(peaks[ch], linestyle=':', color='green')
            if (ch>7):
                axes.axvline(peakWindows[ch-8,0], linestyle=':', color='red')
                axes.axvline(peakWindows[ch-8,1], linestyle=':', color='red')  
        if (Hilbert):
            envelope = getHilbertEnvelopeSingleChannel(deConv_v)
            axes.plot(deConv_t,envelope, alpha=0.7, label="Hilbert Envelope")    
    
#     for ch in range(16):
#         voltage, time = extractChannelWaveform(usefulEvent, ch)
#         if (ch<8):  #VPol case
#             #Perform deconvolution
#             deConv_t,deConv_v = deConvolve_antenna(time, voltage, np.radians(thetaReco[ch]), np.radians(phiReco[ch]), 0)
#         else:  #HPol case
#             #Perform deconvolution
#             deConv_t,deConv_v = deConvolve_antenna(time, voltage, np.radians(thetaReco[ch]), np.radians(phiReco[ch]), 1)
            
#         ax[int(ch/4), ch%4].set_title("Channel " + str(ch))
#         ax[int(ch/4), ch%4].plot(deConv_t,deConv_v, alpha=0.7)
#         if (Hilbert):
#             envelope = getHilbertEnvelopeSingleChannel(deConv_v)
#             ax[int(ch/4), ch%4].plot(deConv_t,envelope, alpha=0.7)
    if (station is not None):
        stationTitle = "Station: " + str(station) + " | "
    else:
        stationTitle = ''
    if (runNumber is not None):
        runTitle = "Run Number: " + str(runNumber) + " | "
    else:
        runTitle = ''
    plt.suptitle(stationTitle + runTitle + "Event Number: " +str(eventNumber) + extraInfo)
   
    if (timeMarkers is not None):
        if (numpy.isscalar(timeMarkers)):
            for ch, axes in zip(channels, ax.flatten()):
                axes.axvline(timeMarkers, color='black', linestyle='--')
        elif (len(timeMarkers)==16):
            for ch, axes in zip(channels, ax.flatten()):
                axes.axvline(timeMarkers[ch], color='black', linestyle='--')    
        # if (numpy.isscalar(timeMarkers)):
        #     for ch in range(16):
        #         ax[int(ch/4), ch%4].axvline(timeMarkers, color='black', linestyle='--')
        # elif (len(timeMarkers)==16):
        #     for ch in range(16):
        #         ax[int(ch/4), ch%4].axvline(timeMarkers[ch], color='black', linestyle='--')
    return fig, ax