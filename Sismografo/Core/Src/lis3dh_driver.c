/*
 * lis3dh_driver.c
 *
 *  Created on: 16/03/2026
 *      Author: joadj
 */

#include "lis3dh_driver.h"

static uint8_t ReadReg(I2C_HandleTypeDef *hi2c, uint8_t reg) {
    uint8_t val = 0;
    HAL_I2C_Mem_Read(hi2c, LIS3DSH_I2C_ADDR, reg, 1, &val, 1, 10);
    return val;
}

static void WriteReg(I2C_HandleTypeDef *hi2c, uint8_t reg, uint8_t val) {
    HAL_I2C_Mem_Write(hi2c, LIS3DSH_I2C_ADDR, reg, 1, &val, 1, 10);
}

uint8_t LIS3DSH_Init(I2C_HandleTypeDef *hi2c) {
    uint8_t who = ReadReg(hi2c, LIS3DSH_WHO_AM_I);
    if (who != 0x3F) return 0;

    // CTRL_REG4: 100Hz ODR, BDU activado, ejes X Y Z habilitados
    WriteReg(hi2c, LIS3DSH_CTRL_REG4, 0x6F);

    // CTRL_REG5: ±2g, ancho de banda 800Hz
    WriteReg(hi2c, LIS3DSH_CTRL_REG5, 0x00);

    return 1;
}

void LIS3DSH_ReadAccel(I2C_HandleTypeDef *hi2c, LIS3DSH_Data *data) {
    uint8_t raw[6];
    // 0x80 activa auto-incremento de registros
    HAL_I2C_Mem_Read(hi2c, LIS3DSH_I2C_ADDR,
                     LIS3DSH_OUT_X_L | 0x80, 1, raw, 6, 10);

    data->x = (int16_t)(raw[0] | (raw[1] << 8));
    data->y = (int16_t)(raw[2] | (raw[3] << 8));
    data->z = (int16_t)(raw[4] | (raw[5] << 8));

    // Sensibilidad: 0.06 mg/LSB en ±2g
    data->ax_g = data->x * 0.00006f;
    data->ay_g = data->y * 0.00006f;
    data->az_g = data->z * 0.00006f;
}
