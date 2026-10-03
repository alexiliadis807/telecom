# import third-party libraries
import numpy as np
import itur

# import local modules
from physics import c
from orbits import Orbit
from isa_simulator import calculate_atmospheric_properties


class Downlink:
    def __init__(
        self,
        electrical_power_available,
        amplifier_efficiency,
        eta_trans_antenna,
        tr_antenna_diameter,
        r_antenna_diameter,
        eta_receiver_antenna,
        selected_band,
        latitude,
        longitude,
        gs_height=0.0,
        availability_percent=99.9,
    ):
        """A method to create a Downlink object

        Args:
            electrical_power_available (_type_): the available electrical power coming from EPS (in Watts)
            amplifier_efficiency (_type_): how much of the given electrical power becomes
            the power for the transmitter antenna.
            eta_trans_antenna: the transmitter antenna efficiency
            tr_antenna_diameter (_type_): the diameter of the transmitter antenna
            r_antenna_diameter: the diameter of the receiver antenna
            eta_receiver_antenna(): the efficiency of the receiver antenna
            selected_band (_type_): at which band is the transmission occuring
            latitude (float): the latitude of the ground station
            longitude (float): the longitude of the ground station
            gs_height (float, optional): the ground station height (in km) from sea level. Default value: 0.0
            (at sea level)
            availability_percent(float): the percentage of the year that the satellite link
            is operating within acceptable performance thresholds.
        """
        # the bandwidths are found based on values reported by the NASA
        self.bands = {
            "S-band": {"f_central": 3e9, "bandwidth": 2.0e6},
            "X-band": {"f_central": 10e9, "bandwidth": 10.0e6},
            "Ka-band": {"f_central": 33.5e9, "bandwidth": 100.0e6},
        }

        self.selected_band = selected_band
        self.f = self.bands[selected_band]["f_central"]

        self.electrical_power_available = electrical_power_available
        self.amplifier_efficiency = amplifier_efficiency
        self.transmission_power = electrical_power_available * amplifier_efficiency
        self.eta_trans_antenna = eta_trans_antenna

        self.latitude = latitude
        self.longitude = longitude
        self.gs_height = gs_height

        self.availability_percent = availability_percent
        self.p_exceedance = 100.0 - self.availability_percent

        self.tr_antenna_diameter = tr_antenna_diameter
        self.r_antenna_diameter = r_antenna_diameter
        self.eta_receiver = eta_receiver_antenna

    def choose_band(self, band):
        if band not in self.bands:
            message = "Please select another band (S-band, X-band, Ka-band)"
            return message
        else:
            self.selected_band = band
            self.f_central = self.bands[band]["f_central"]
            self.bandwidth = self.bands[band]["bandwidth"]
            self.wavelength = c / self.f_central
            return f"Successfully switched to {band}"

    def determine_free_space_path_loss(self, distance: float):
        """a method to calculate the FSPL.

        Args:
            distance (float): the distance (in m) to the telescope/satellite.
        Returns:
            float: the fspl in dB.
        """
        fspl = ((4 * np.pi * distance * self.f_central) / (c)) ** 2

        # convert to dB
        self.fspl_dB = 10 * np.log10(fspl)

        return self.fspl_dB

    def determine_gas_attenuation(self, elevation_angle_deg):
        """
        we will use the modules from itur to compute the attenuations.
        the following method will give us the gas attenuation based on the lower limit we set
        this will be set at h = 0 km based on standard atmosphere values (as given by ITU-R P.676).
        the gas attenuation is an integral of γ(h) (how many dB are lost for each km) along the slant path. The model applies the method
        of equivalent heights to calculate precise amount of loss in dB along the diagonal path of our
        satellite
        """
        # first calculate local atmospheric properties based on coords and height
        h_m = self.gs_height * 1000
        T_k, P_Pa, rho_isa = calculate_atmospheric_properties(h_m)
        P_hPa = P_Pa / 100.0

        # now from  ITU-R  P.835-7 we also know
        rho_sea = 7.5  # g/m^3 (density of H2O(g) at sea level)
        rho_exp = rho_sea * np.exp(
            -self.gs_height / 2
        )  # where 2 is the scale height, H_0, of H2O(g) in atmosphere (this is 2km)

        # lowest limit based on mixing ratio (equation 8 from "ITU-R P.835-7", pg.5):
        rho_min = 2e-6 * ((P_hPa * 216.7) / T_k)

        # finally find the density of the water vapor
        rho_wv = max(rho_exp, rho_min)

        # and calculate the gas attenuation
        gas_attenuation = itur.gaseous_attenuation_slant_path(
            self.f_central / 1e9, elevation_angle_deg, rho=rho_wv, P=P_hPa, T=T_k
        )

        res_array = np.asarray(gas_attenuation)
        return float(np.squeeze(res_array))

    def determine_rain_attenuation(self, elevation_angle_deg):
        """Computes Rain Attenuation (A_rain) according to ITU-R P.838 models.

        Args:
            elevation_angle_deg (float): The angle (in degrees) the satellite is wrt the local horizon.
            availability_percent (float, optional): Desired link availability. Defaults to the
            standard and sweet spot of 99.9. This percentage would lead to 0.1% exceedance probability.
        """
        # calculate rain attenuation from itur
        rain_attenuation = itur.rain_attenuation(
            lat=self.latitude,
            lon=self.longitude,
            f=self.f_central / 1e9,
            el=elevation_angle_deg,
            p=self.p_exceedance,
            hs=self.gs_height,
        )

        res_array = np.asarray(rain_attenuation)
        return float(np.squeeze(res_array))

    def determine_cloud_attenuation(self, elevation_angle_deg):
        cloud_attenuation = itur.cloud_attenuation(
            lat=self.latitude,
            lon=self.longitude,
            el=elevation_angle_deg,
            f=self.f_central / 1e9,
            p=self.p_exceedance,
        )
        res_array = np.asarray(cloud_attenuation)
        return float(np.squeeze(res_array))

    def determine_scintillation_attenuation(self, elevation_angle_deg):
        """This method calculates the Tropospheric Scintillation Loss (A_scint) in db
        following the method as presented in ITU-R P.618 (section 2.4)
        """
        scint_attenuation = itur.scintillation_attenuation(
            lat=self.latitude,
            lon=self.longitude,
            f=self.f_central / 1e9,
            el=elevation_angle_deg,
            p=self.p_exceedance,
            D=self.r_antenna_diameter,
            eta=self.eta_receiver,
        )
        res_array = np.asarray(scint_attenuation)
        return float(np.squeeze(res_array))

    def determine_total_atmospheric_loss(self, elevation_angle_deg):
        A_gas = self.determine_gas_attenuation(elevation_angle_deg)
        A_rain = self.determine_rain_attenuation(elevation_angle_deg)
        A_cloud = self.determine_cloud_attenuation(elevation_angle_deg)
        A_scint = self.determine_scintillation_attenuation(elevation_angle_deg)

        # based on ITU-R P.618 (section 2.5)
        A_total = float(A_gas + np.sqrt((A_rain + A_cloud) ** 2 + A_scint**2))

        attenuations = {
            "A_gas": A_gas,
            "A_rain": A_rain,
            "A_cloud": A_cloud,
            "A_scint": A_scint,
            "A_total": A_total,
        }
        return attenuations

    def calculate_antenna_gain(self, eta: float, diameter: float) -> float:
        """A method to calculate the gain of the transmitter and the receiver antenna.

        Args:
            eta (float): the efficiency of the antenna
            diameter (float): the diameter of the antenna

        Returns:
            gain(float): the gain of the antenna in dBi (decibels relative to isotropic)
        """
        gain = eta * (np.pi * diameter / self.wavelength) ** 2
        gain_dB = 10 * np.log10(gain)
        return gain_dB

    def calculate_eirp(self, cable_losses=0.5):
        """A method to calculate the effective isotropically radiated power (eirp) by the
        transmitter antenna.

        Args:
            cable_losses (float, optional): The losses inside the cables that connect
            the amplifier to the antenna (those losses may include heat dissipation due to
            ohmic resistance, energy absorbed by the dielectrices, etc.). Defaults to 0.5 dB.
        Returns:
            eirp_dBW: the eirp in dBW
        """

        transmission_power_dBW = 10 * np.log10(self.transmission_power)
        transmission_gain = self.calculate_antenna_gain(
            self.eta_trans_antenna, self.tr_antenna_diameter
        )
        eirp_dBW = transmission_power_dBW + transmission_gain - cable_losses

        return eirp_dBW

    def calculate_g_over_t(
        self, elevation_angle_deg: float, T_mr: float = 275.0, NF_dB: float = 1.5
    ):
        """A method to calculate the ratio of the gain of the receiving antenna to the noise
        that covers the signal. This will be used as a measure of the efficiency of the
        ground station.

        Args:
            elevation_angle_deg (float): the elevation angle of the satellite with respect to the
            local horizon (in degrees).
            T_mr (float, optional): The atmospheric mean radiating temperature as
            specified by ITU-R P.618-14 (in section: "Noise Temperature"). Defaults to 275.
            NF_dB (float, optional): The noise figure in decibels, which acts as a figure
            of merit. Defaults to 1.5.
        """
        # STEP 0: Employ the methods from the Downlink class to determine
        # a) The total atmospheric attenuation excluding scintillation fading
        # b) The gain of the receiving antenna
        attenuations = self.determine_total_atmospheric_loss(elevation_angle_deg)
        A_total = attenuations["A_total"]
        A_scint = attenuations["A_scint"]
        A = A_total - A_scint

        reception_gain = self.calculate_antenna_gain(
            self.eta_receiver, self.r_antenna_diameter
        )

        # STEP 1: Use equation 69 from ITU-R P.618-14 to calculate Tsky
        T_sky = T_mr * (1 - 10 ** (-A / 10)) + 2.7 * (10 ** (-A / 10))

        # STEP 2: Calculate T_rx
        F_n = 10 ** (NF_dB / 10)
        T_rx = 290 * (F_n - 1)

        # STEP 3: Calculate T_sys:
        T_sys = T_sky + T_rx

        # STEP 4: Calculate G/T
        g_over_t = reception_gain - 10 * np.log10(T_sys)

        return g_over_t

    def calculate_c_over_n0(
        self, distance: float, elevation_angle_deg: float, other_losses_dB: float = 1.5
    ) -> float:
        """A method to calculate the Carrier-to-Noise density ratio in dB-Hz.

        Args:
            distance (float): Slant Range Distance to the satellite in meters.
            elevation_angle_deg (float): the elevation angle of the satellite with respect to the
            local horizon (in degrees).
            other_losses_dB (float, optional): Pointing, polarization and hardware mismatch losses in dB.
            Defaults to 1.5.

        Returns:
            float: Carrier-to-Noise density ratio in dB-Hz.
        """
        # STEP 0: Employ the methods from the Downlink class
        eirp = self.calculate_eirp()
        fspl = self.determine_free_space_path_loss(distance)
        atmospheric_attenuation = self.determine_total_atmospheric_loss(
            elevation_angle_deg
        )["A_total"]
        g_over_t = self.calculate_g_over_t(elevation_angle_deg)

        # boltzmann constant in dBW/(K*Hz)
        k_dB = -228.6

        # STEP 1: Calculate c_over_n0 (in line with "Satellite Communications" by T. Pratt and J. Allnutt (chapter 4)):
        cn0 = eirp - fspl - atmospheric_attenuation - other_losses_dB + g_over_t - k_dB
        self.cn0 = cn0

        return cn0

    def determine_data_rate(
        self, daily_data_vol: float, overhead_data: float, t_pass_hours: int = 4
    ):
        """A method to calculate the data rate required for the downlink (including a margin for
        the overhead data)

        Args:
            daily_data_vol (float): the data volume to be transmitted per day (in Gbit)
            overhead_data (float): makes sure data received by the ground station can be clearly
            read (includes FEC, CCSDS space packets)
            t_pass_hours (int, optional): how many hours the ground contact last.
                        Defaults to 4 hours.
        Returns:
            float: the data rate required (accounting for the margin (overhead data)).
        """
        t_pass_sec = t_pass_hours * 3600
        R_raw = (daily_data_vol * (10**9)) / t_pass_sec
        R_req = R_raw * (1 + overhead_data)
        self.data_rate_required = R_req

        return R_req

    def calculate_eb_over_n0(
        self,
        distance: float,
        elevation_angle_deg: float,
        daily_data_vol: float,
        overhead_data: float,
        t_pass_hours: int = 4,
    ) -> float:
        """Calculates achieved Eb/N0 (in dB) directly from inputs."""
        cn0 = self.calculate_c_over_n0(distance, elevation_angle_deg)
        r_req = self.determine_data_rate(daily_data_vol, overhead_data, t_pass_hours)

        R_dB = 10 * np.log10(r_req)
        eb_over_n0 = cn0 - R_dB
        return eb_over_n0

    def calculate_shannon_capacity(
        self, distance: float, elevation_angle_deg: float
    ) -> float:
        """Calculates theoretical Shannon Channel Capacity (in bps)."""
        cn0 = self.calculate_c_over_n0(distance, elevation_angle_deg)

        SNR_dB = cn0 - 10 * np.log10(self.bandwidth)
        SNR_lin = 10 ** (SNR_dB / 10)
        shannon_capacity = self.bandwidth * np.log2(1 + SNR_lin)

        return shannon_capacity

    def calculate_link_margin(self, distance:float, elevation_angle_deg: float, 
                              daily_data_vol: float, overhead_data:float, 
                              eb_n0_threshold_dB:float=3.0, t_pass_hours:int = 4):
        """Calculates the Link Margin

        Args:
            distance (float): Slant range distance to the satellite in meters.
            elevation_angle_deg (float): Elevation angle in degrees.
            daily_data_vol (float): Daily data volume to be transmitted (in Gbit).
            overhead_data (float): Makes sure data received by the ground station can be clearly
            read (includes FEC, CCSDS space packets).
            eb_n0_threshold_dB (float, optional): According to ECSS-E-ST-50-05C (section 4.2) "The EIRP 
            transmitted from the Earth station shall be selected to allow for 
            a margin of 3 dB on the link budget". Defaults to 3.0.
            t_pass_hours (int, optional): How many hours the ground contact last. Defaults to 4.

        Returns:
            float: Link Margin in dB.
            str: a comment on whether the link budget is compliant.
        """

        eb_n0_achieved = self.calculate_eb_over_n0(distance, elevation_angle_deg, daily_data_vol, 
                                                   overhead_data, t_pass_hours)

        link_margin = eb_n0_achieved - eb_n0_threshold_dB
        
        self.required_margin = 3.0 
        if link_margin >= self.required_margin:
            statement = "STATUS: PASSED (Link Budget Compliant)"
        else:
            statement = f"STATUS: FAILED (Deficit of {round(self.required_margin - np.abs(link_margin),2)} dB)"

        return link_margin, statement
    
        

