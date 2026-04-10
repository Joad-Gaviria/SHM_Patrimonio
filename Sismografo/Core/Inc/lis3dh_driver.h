/*
 * lis3dh_driver.h
 *
 *  Created on: 16/03/2026
 *      Author: joadj
 */

#ifndef INC_LIS3DH_DRIVER_H_
#define INC_LIS3DH_DRIVER_H_

#include "stm32h7xx_hal.h"

// Dirección I2C correcta para LIS3DSH
#define LIS3DSH_I2C_ADDR    (0x1D << 1)

// Registros del LIS3DSH (distintos al LIS3DH)
#define LIS3DSH_WHO_AM_I    0x0F  // Devuelve 0x3F
#define LIS3DSH_CTRL_REG4   0x20  // ODR, ejes, BDU
#define LIS3DSH_CTRL_REG5   0x24  // Escala, ancho de banda
#define LIS3DSH_OUT_X_L     0x28

typedef struct {
    int16_t x, y, z;
    float ax_g, ay_g, az_g;
} LIS3DSH_Data;

uint8_t LIS3DSH_Init(I2C_HandleTypeDef *hi2c);
void    LIS3DSH_ReadAccel(I2C_HandleTypeDef *hi2c, LIS3DSH_Data *data);

#endif /* INC_LIS3DH_DRIVER_H_ */
