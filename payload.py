'''
A module to insert the data for the payload
'''

class Payload:
    def __init__(self,
                 swath_width_angle,
                 pixel_size,
                 bits_per_pixel, duty_cycle, 
                 downlink_time):

        '''
        swath_width_angle: the total ground width the sensor covers in a single pass
        pixel_size: the angular resolution represented by a single detector pixel 
        bits_per_pixel: the number of digital data bits used to represent the signal intensity
        of each pixel
        duty_cycle: the fraction of time the payload active instrument is powered on and is 
        taking data relative to the total orbit duration
        downlink_time: time available to transmit payload data to the ground station
        '''
        
        self.swath_width_angle = swath_width_angle
        self.pixel_size = pixel_size
        self.bits_per_pixel = bits_per_pixel
        self.duty_cycle = duty_cycle
        self.downlink_time = downlink_time