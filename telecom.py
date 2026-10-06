# import third-party libraries
import numpy as np
import itur

# import local modules
from payload import Payload


class Downlink:
    def __init__(
        self,
        f_downlink,
        sc_tx_power,
        gs_tx_power,
        sc_antenna_diameter,
        gs_antenna_diameter,
        turn_around_ratio,
        orbit_altitude,
        elongation_angle,
        pointing_offset_angle,
        required_uplink_dr,
        payload: Payload,
        required_bit_error_rate,
        elevation_angle=10,
        loss_factor_tx=0.8,
        loss_factor_rx=0.7,
    ):
        """A method to create a Downlink object

        Args:
            sc_tx_power (int): the transmitter power of the s/c.
            gs_tx_power (int): the transmitter power of the ground station.
            loss_factor_tx (float): the transmitter antenna loss factor. It defaults to 0.8
            f_downlink(float): the downlink frequency (in GHz)
            sc_antenna_diameter (float): the diameter of the antenna of the s/c.
            gs_antenna_diameter (float): the diameter of the antenna of the ground station.
            loss_factor_rx (float): the receiver antenna loss factor. It defaults to 0.7
            turn_around_ratio (float): the uplink to downlink frequency ratio
            orbit_altitude(int): the altitude of the s/c in km.
            elongation_angle(int): angle between s/c - Sun line and Earth - Sun line in deg.
            pointing_offset_angle(float): the angle between the actual direction the s/c is
            pointing and the desired direction of the target.
            required_uplink_dr(int): the data rate required for uplink
            required_bit_error_rate(int): the maximum acceptable ratio of erroneous
            bits to total transmitted bits.
            payload (Payload): Instance of Payload containing spacecraft payload specification attributes.
            elevation_angle(int): The min elevation angle wrt local horizon. Assumed equal to
            10degrees for default.
        """

        # Store the payload object
        self.payload = payload

        self.tx_loss = loss_factor_tx
        self.rx_loss = loss_factor_rx

        self.f_downlink = f_downlink

        self.sc_tx_power = sc_tx_power
        self.gs_tx_power = gs_tx_power
        self.turn_around_ratio = turn_around_ratio

        self.sc_antenna_diameter = sc_antenna_diameter
        self.gs_antenna_diameter = gs_antenna_diameter

        # assume both antennas to be parabolic
        self.antenna_efficiency = 0.55

        self.orbital_alt = orbit_altitude
        self.elongation_angle = elongation_angle
        self.elevation_angle = elevation_angle

        # based on this the pointing loss during downlink will be estimated.
        self.pointing_offset_e = pointing_offset_angle

        self.uplink_dr_required = required_uplink_dr

        self.bit_error_rate = required_bit_error_rate

        # finally define the constants to be used
        self.c = 3 * 10**8  # m/s
        self.earth_radius = 6371.0  # km
        self.AU_to_km = 149597870.7  # km/AU
        self.moon_earth_distance = 384400  # km

    def get_rf_params(self, mode: str):
        if mode.lower() == "downlink":
            return (
                self.f_downlink,
                self.sc_tx_power,
                self.sc_antenna_diameter,
                self.gs_antenna_diameter,
            )
        else:
            f_uplink = self.f_downlink * self.turn_around_ratio
            return (
                f_uplink,
                self.gs_tx_power,
                self.gs_antenna_diameter,
                self.sc_antenna_diameter,
            )

    def determine_pointing_loss(self, mode):
        """
        Args:
            mode: the type of communication mode (uplink, downlink).
        Returns:
            L_pointing: the pointing loss in [dB]
        """
        f, _, _, _ = self.get_rf_params(mode)
        # ASSUMPTION: PARABOLIC ANTENNA (as this is the only type of antenna considered
        # in the study cases)
        d_sc = self.sc_antenna_diameter
        half_power_angle = 21 / (f * d_sc)

        # calculate the pointing loss in dB
        L_pointing = 12 * (self.pointing_offset_e / half_power_angle) ** 2

        return L_pointing

    def determine_free_space_path_loss(self, mode, type, d_S):
        """a method to calculate the FSPL.

        Args:
            mode: the type of communication mode (uplink, downlink).
            type(str): defines the type of mission (lunar orbit, earth orbit, interplanetary mission)
            d_S (float): defines the distance between the s/c and the Sun and may be
            considered equal to the distance of the Planet (which the s/c orbit) to the Sun. It is in km.
        Returns:
            float: the fspl in dB.
        """
        f, tx_power, d_tx, d_rx = self.get_rf_params(mode)
        if type == "Earth Orbit":
            R_e = self.earth_radius
            h = self.orbital_alt
            alpha = np.radians(self.elevation_angle)
            S = R_e * (
                np.sqrt((((h + R_e) / R_e) ** 2) - (np.cos(alpha) ** 2)) - np.sin(alpha)
            )

        elif type == "Lunar Orbit":
            # Assume distance between Earth and Moon is much bigger than the orbital attitude.
            S = self.moon_earth_distance

        else:
            d_ES = self.AU_to_km #km
            theta_ES = np.radians(self.elongation_angle)
            S = np.sqrt((d_ES**2) + (d_S**2) - 2 * d_ES * d_S * np.cos(theta_ES))

        S_m = 1000 * S
        fspl = ((4 * np.pi * S_m * f * (10**9)) / (self.c)) ** 2

        # convert to dB
        self.fspl_dB = 10 * np.log10(fspl)

        return self.fspl_dB

    def determine_atmospheric_attenuation(self, mode: str):
        """Calculates atmospheric attenuation (in dB).
        Args:
            mode (str): Communication mode ('uplink' or 'downlink').

        Returns:
            float: Atmospheric attenuation in dB.
        """
        # Since all links connect to an Earth Ground Station,

        # Earth slant-path attenuation always applies.
        T_k = 288.15  # K (Standard ISA Sea Level)
        P_hPa = 1013.25  # hPa
        rho_sea = 7.5  # g/m^3

        f, _, _, _ = self.get_rf_params(mode)
        elevation = max(self.elevation_angle, 5.0)  # Bound for ITU-R model validity

        gas_attenuation = itur.gaseous_attenuation_slant_path(
            f, elevation, rho=rho_sea, P=P_hPa, T=T_k
        )
        return float(np.squeeze(np.asarray(gas_attenuation)))

    def calculate_transmission_gain(self, mode):
        """A method to calculate the gain of the transmitter and the receiver antenna.

        Returns:
            gain(float): the gain of the antenna in dBi (decibels relative to isotropic)
        """
        f, _, d_tx, _ = self.get_rf_params(mode)
        gain = self.antenna_efficiency * ((np.pi * d_tx * f * (10**9)) / self.c) ** 2
        gain_dB = 10 * np.log10(gain)
        return gain_dB

    def calculate_reception_gain(self, mode: str) -> float:
        """Calculates receiving antenna gain in dBi."""
        f, _, _, d_rx = self.get_rf_params(mode)
        gain = self.antenna_efficiency * ((np.pi * d_rx * (f * 1e9)) / self.c) ** 2
        return float(10 * np.log10(gain))

    def calculate_eirp(self, mode):
        """A method to calculate the effective isotropically radiated power (eirp) by the
        transmitter antenna.

        Args:
            mode: set the mode to either uplink or downlink, so that you get the correct EIRP.
        Returns:
            eirp_dBW: the eirp in dBW
        """
        f, tx_power, d_tx, d_rx = self.get_rf_params(mode)

        transmission_power_dBW = 10 * np.log10(tx_power)

        # convert tx_loss in dB
        tx_loss_dB = -10 * np.log10(
            self.tx_loss
        )  # self.tx_loss < 1 and thus -10*log10(self.tx_loss) = +0.969 dB
        # AN ASSUMPTION MADE HERE IS THAT THE EFFICIENCY FOR TRANSMITTER AND RECEIVER
        # ANTENNAS ARE THE SAME FOR BOTH DOWNLINK AND UPLINK (i.e. independent of mode)
        transmission_gain = self.calculate_transmission_gain(mode)
        eirp_dBW = transmission_power_dBW + transmission_gain - tx_loss_dB

        return eirp_dBW

    def calculate_g_over_t(self, mode, type, T_mr: float = 275.0, NF_dB: float = 1.5):
        """A method to calculate the ratio of the gain of the receiving antenna to the noise
        that covers the signal. This will be used as a measure of the efficiency of the
        ground station.

        Args:
            T_mr (float, optional): The atmospheric mean radiating temperature as
            specified by ITU-R P.618-14 (in section: "Noise Temperature"). Defaults to 275K,
            in the absence of local data.
            NF_dB (float, optional): The noise figure in decibels, which acts as a figure
            of merit. It is assumed equal to 1.5.
            mode(str): set the mode to either uplink or downlink
            type(str): type of mission (e.g. "Earth Orbit", "Lunar Orbit", etc.)
        """
        # define the reception gain
        reception_gain = self.calculate_reception_gain(mode)  # dB

        if mode.lower() == "downlink":
            # we assume clear-sky as no weather conditions or ground station coordinates are provided.
            # so A = Ag (gaseous absorption), in the ITU-R formula.
            A = self.determine_atmospheric_attenuation(mode)

            """STEP 1: Use equation 69 from ITU-R P.618-14 to calculate Tant. In fact the equation 
            calculates Tsky (sky noise temperature at the ground station antenna), but we assume it equal to Tant.
            The effect of this assumption is that it neglects ground brightness 
            temperature picked up by antenna sidelobes."""
            T_ant = T_mr * (1 - 10 ** (-A / 10)) + 2.7 * (10 ** (-A / 10))

        else:
            # Spacecraft points at Earth / target planet -> Planet background noise (assume Tant equal to Earth
            # brightness temperature)
            T_ant = 290.0

        # STEP 2: Calculate T_rx
        F_n = 10 ** (NF_dB / 10)
        T_n = 290 * (F_n - 1)
        T_cable = 290 * (1 - self.rx_loss) / self.rx_loss

        # STEP 3: Calculate T_sys.:
        T_sys = T_ant + T_n + T_cable

        # STEP 4: Calculate G/T
        g_over_t = reception_gain - 10 * np.log10(T_sys)

        return g_over_t

    def calculate_c_over_n0(self, mode, type, d_S, other_losses_dB=1.5):
        """A method to calculate the Carrier-to-Noise density ratio in dB-Hz.

        Args:
            mode(str): set the mode to either uplink or downlink
            elevation_angle_deg (float): the elevation angle of the satellite with respect to the
            local horizon (in degrees).
            type(str): defines the type of mission (lunar orbit, earth orbit, interplanetary mission)
            d_S (float): defines the distance between the s/c and the Sun in m and may be
            considered equal to the distance of the Planet (which the s/c orbit) to the Sun.
            other_losses_dB (float, optional): Polarization and hardware mismatch losses in dB.
            Defaults to 1.5.

        Returns:
            float: Carrier-to-Noise density ratio in dB-Hz.
        """
        # STEP 1: Employ the methods from the Downlink class
        eirp = self.calculate_eirp(mode)
        fspl = self.determine_free_space_path_loss(mode, type, d_S)
        pointing_loss = self.determine_pointing_loss(mode)
        atmospheric_attenuation = self.determine_atmospheric_attenuation(mode)

        g_over_t = self.calculate_g_over_t(mode, type)

        # boltzmann constant in dBW/(K*Hz)
        k_dB = -228.6

        # STEP 2: Calculate c_over_n0 (in line with "Satellite Communications" by T. Pratt and J. Allnutt (chapter 4)):
        cn0 = (
            eirp
            - fspl
            - atmospheric_attenuation
            - pointing_loss
            - other_losses_dB
            + g_over_t
            - k_dB
        )
        self.cn0 = cn0

        return cn0

    def determine_data_rate(self, mu, body_radius, mode):
        """
        A method to calculate the payload data rate.
        Args:
            mu(float): the gravitational parameter of body that is being orbited (in m^3/s^2).
            body_radius(float): the radius of the body (e.g. Mars, Mercury, Earth) being orbitted in km.
        """
        if mode.lower() == "downlink":
            # Duty cycle
            D_c = self.payload.duty_cycle

            # Ratio of actual downlink time to total available tim
            T_DL = self.payload.downlink_time / 24.0

            # bits per pixel:
            B_p = self.payload.bits_per_pixel

            # swath width:
            S_w_deg = self.payload.swath_width_angle
            S_w_rad = np.radians(S_w_deg)

            # pixel size:
            pixel_size_arcmin = self.payload.pixel_size
            pixel_size_rad = (pixel_size_arcmin / 60) * (np.pi / 180)

            # pixels per line:
            N_p = S_w_rad / pixel_size_rad

            # Ground pixel size in m:
            P_s = pixel_size_rad * (1000 * self.orbital_alt)

            # Calculate ground speed assuming circular orbit
            V_orbit = np.sqrt(mu / (1000 * (body_radius + self.orbital_alt)))  # m/s
            V_ground = V_orbit * (
                (body_radius) / (body_radius + self.orbital_alt)
            )  # m/s

            # Calculate lines per second:
            f_line = V_ground / P_s

            # calculate the generated data rate
            R_G = B_p * N_p * f_line

            # Now calculate the required data rate
            R = R_G * (D_c / T_DL)

        else:
            R = float(self.uplink_dr_required)

        return R

    def calculate_eb_over_n0(self, mode, type, d_S, body_radius, mu) -> float:
        """Calculates achieved Eb/N0 (in dB) directly from inputs.
        The arguments are the same as the ones in the other methods.
        """

        cn0 = self.calculate_c_over_n0(mode, type, d_S)
        r_req = self.determine_data_rate(mu, body_radius, mode)

        # covert data rate to dB
        R_dB = 10 * np.log10(r_req)

        # calculate eb_over_n0 in dB
        eb_over_n0 = cn0 - R_dB

        return eb_over_n0

    def calculate_link_margin(
        self, mode, type, d_S, mu, body_radius, eb_n0_threshold_dB=3.0
    ):
        """Calculates the Link Margin

        Args:
            mode(str): set the mode to either uplink or downlink
            d_S (float): defines the distance between the s/c and the Sun in m and may be
            considered equal to the distance of the Planet (which the s/c orbit) to the Sun.
            type(str): defines the type of mission (lunar orbit, earth orbit, interplanetary mission)
            eb_n0_threshold_dB (float, optional): According to ECSS-E-ST-50-05C (section 4.2) "The EIRP
            transmitted from the Earth station shall be selected to allow for
            a margin of 3 dB on the link budget". Defaults to 3.0.

        Returns:
            float: Link Margin in dB.
            str: a comment on whether the link budget is compliant.
        """

        eb_n0_achieved = self.calculate_eb_over_n0(mode, type, d_S, body_radius, mu)

        link_margin = eb_n0_achieved - eb_n0_threshold_dB

        if link_margin >= 0.0:
            statement = "STATUS: PASSED (Link Budget can be closed)"
        else:
            statement = (
                f"STATUS: FAILED (Deficit of {round(abs(link_margin), 2)} dB). "
                "Consider changing the transmitter power or using a more powerful antenna."
            )

        return link_margin, statement
