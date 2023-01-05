import sys

import xlrd
from PySide2.QtUiTools import QUiLoader
from PySide2.QtWidgets import QApplication, QMessageBox, QColorDialog, QWidget, QMainWindow
from PySide2.QtCore import QFile, QCoreApplication, Qt, QSize, QObject, QEvent
from PySide2.QtGui import QIcon, QMovie
import PySide2.QtXml

import numpy as np
import math
from scipy.fft import fft
import matplotlib.pyplot as plt
from time import sleep

import os


class CoastalDataProc(QWidget):
    time_step = 0.01
    WaveformColor = 'red'
    DrawPointsColor = 'red'
    MaxHeightPointColor = 'green'
    MaxHeightLineColor = 'black'
    MinHeightPointColor = 'blue'
    MinHeightLineColor = 'black'
    ZeroPointsColor = 'blue'
    datapointsColor = 'red'
    bShowColorSetting = False

    sy1_H = 0.3  # 实验一水深0.3米

    def __init__(self):
        super(CoastalDataProc, self).__init__()
        self.ui = QUiLoader().load("GUI.ui")

        self.load_sy1_files()
        self.load_sy2_files()
        self.load_sy3_files()

        self.ui.WaveformColor_btn.setStyleSheet('background: {};'.format('#ff0000'))
        self.ui.DrawPointsColor_btn.setStyleSheet('background: {};'.format('#ff0000'))
        self.ui.MaxHeightPointColor_btn.setStyleSheet('background: {};'.format('#00ff00'))
        # self.ui.MaxHeightLineColor_btn.setStyleSheet('background: {};'.format('#000000'))
        self.ui.MinHeightPointColor_btn.setStyleSheet('background: {};'.format('#0000ff'))
        # self.ui.MinHeightLineColor_btn.setStyleSheet('background: {};'.format('#000000'))
        self.ui.ZeroPointsColor_btn.setStyleSheet('background: {};'.format('#0000ff'))
        # self.ui.datapointsColor_btn.setStyleSheet('background: {};'.format('#ff0000'))

        self.gif = QMovie('res\\duck.gif')
        self.ui.label_gif.setMovie(self.gif)
        self.gif.setScaledSize(QSize(32, 32))

        self.gif.start()
        self.ui.show_colorsetting.setIcon(QIcon("res\\color.jpg"))
        self.ui.groupBox_4.hide()
        # a=colorchooser.askcolor()
        # print(a)

    def close_Event(self):
        plt.close('all')  # 程序退出时关闭所有窗口

    def load_sy1_files(self):
        for root, dirs, files in os.walk("datas\第1次海动实验"):
            if files:
                n_files = []
                f_name = root.split('\\')[-2:]
                n_files.append('\\'.join(f_name))
                self.ui.filelist_sy1.addItems(n_files)

    def load_sy2_files(self):
        for root, dirs, files in os.walk("datas\第2次海动实验"):
            if files:
                n_files = []
                f_name = root.split('\\')[-1]
                for i in files:
                    n_files.append(f_name + '\\' + i)
                self.ui.filelist_sy2.addItems(n_files)

    def load_sy3_files(self):
        sy3_filelist = []
        for root, dirs, files in os.walk("datas\第3次海动实验"):
            if files:
                if '波高测量' in root:
                    n_files = []
                    f_name = root.split('\\')[-3:-1]
                    n_files.append('\\'.join(f_name))
                    self.ui.filelist_sy3.addItems(n_files)

    def btn_StartProc_sy1_click(self):
        reply = QMessageBox.question(self, '提醒',
                                     "请确认造波仪输入波高已经正确填写：\nH1 = {}cm\nH2 = {}cm\nH3 = {}cm".format(
                                         self.ui.sy1_contrastform_H_1.value(), self.ui.sy1_contrastform_H_2.value(),
                                         self.ui.sy1_contrastform_H_3.value()), QMessageBox.Yes | QMessageBox.No,
                                     QMessageBox.No)
        if reply == QMessageBox.No:
            return
        plt.close('all')
        self.ui.sy1_WaveFormOut.setText('')
        for root, dirs, files in os.walk("datas\第1次海动实验\\" + self.ui.filelist_sy1.currentText() + '\\'):
            if files:
                self.sy1_rawdata = []
                for sig_file in files:
                    with open(root + sig_file, encoding='utf-8') as file:
                        sy1rawdata = [[adata.strip() for adata in data.split('    ')] for data in
                                      [data.strip() for data in file.readlines()]]
                    del (sy1rawdata[0:2])
                    del (sy1rawdata[-2:])
                    self.sy1_rawdata.append(sy1rawdata)
        self.sy1_coldatas = []
        for i_data in self.sy1_rawdata:
            sy1coldata = [[float(data[i]) for data in i_data] for i in range(len(i_data[0]))]
            self.sy1_coldatas.append(sy1coldata)
        self.sy1_timeline = list(self.floatRange(0, len(i_data), 1, 2))

        self.sy1_smoothdatas = []
        for sy1_t in self.sy1_coldatas:
            sy1_smooth_t = []
            for sy1_j in sy1_t:
                sy1_smooth_t.append(self.moving_average(sy1_j, self.ui.spinBox_smoothWindowSize_sy1.value()))
            self.sy1_smoothdatas.append(sy1_smooth_t)
        drawNumber = self.ui.spinBox_number_sy1.value()
        drawIndex = self.ui.spinBox_index_sy1.value() - 1

        for i in range(drawIndex, min(drawNumber + drawIndex, len(self.sy1_coldatas))):
            fig = plt.figure()
            fig.canvas.set_window_title("<规则波实验> 第{}次实验".format(i + 1))
            ax_waveform = plt.subplot2grid((2, 1), (0, 0), colspan=1, rowspan=1)
            self.sy1_drawWaveform(i, ax_waveform,self.ui.checkBox_sy1_bSmoothData.isChecked(),self.ui.checkBox_sy1_bDrawPoints.isChecked())

            ax_waveform = plt.subplot2grid((2, math.floor(self.ui.sy1_contrastform_Position.value())), (1, 0), colspan=1, rowspan=1)
            self.sy1_ContrastWaveForm(i, ax_waveform)

        plt.rcParams['font.sans-serif'] = ['SimHei']
        plt.rcParams['axes.unicode_minus'] = False
        plt.show()

    def btn_StartProc_sy2_click(self):
        plt.close('all')
        # msg = QMessageBox()
        # msg.setText('数据:{}\n'.format(self.ui.filelist_sy2.currentText()))
        # msg.exec_()
        with open("datas\第2次海动实验\\" + self.ui.filelist_sy2.currentText(), encoding='utf-8') as file:
            self.content = [[adata.strip() for adata in data.split('    ')] for data in
                            [data.strip() for data in file.readlines()]]
        del (self.content[0:2])
        del (self.content[-2:])

        self.coldata = [[float(data[i]) for data in self.content] for i in range(len(self.content[0]))]
        self.timeline = list(self.floatRange(0, len(self.content), 1, 2))

        self.smooth_data = []
        for smoothIndex in range(len(self.coldata)):
            self.smooth_data.append(
                self.moving_average(self.coldata[smoothIndex], self.ui.spinBox_smoothWindowSize.value()))  # 滑动平均降噪

        drawIndex = self.ui.spinBox_index.value() - 1
        selectDrawNumber = min(self.ui.spinBox_number.value(), len(self.coldata))
        drawNumber = self.ui.spinBox_number.value()
        subWaveHeightList = []  # 储存每个传感器的每个周期的波高
        for i in range(drawIndex, min(drawIndex + drawNumber, len(self.coldata))):
            fig = plt.figure()
            fig.canvas.set_window_title("<随即波实验>波高传感器#{}".format(i + 1))
            ax1 = plt.subplot2grid((2, 7), (0, 0), colspan=7, rowspan=1)
            subWaveHeightList.append(self.drawWaveform(i, ax1, self.ui.checkBox_bRawData.isChecked(),
                                                       self.ui.checkBox_bDrawPoints.isChecked(),
                                                       self.ui.checkBox_bDrawMaxHeight.isChecked(),
                                                       self.ui.checkBox_bDrawMinHeight.isChecked(),
                                                       self.ui.checkBox_bDrawZeroPoints.isChecked()))

            canv_WaveHeight = plt.subplot2grid((2, 7), (1, 0), colspan=2, rowspan=1)
            aH = self.drawWaveHeightStatistics(i, canv_WaveHeight, subWaveHeightList[i - drawIndex])

            canv_RayleighDistribution = plt.subplot2grid((2, 7), (1, 4), colspan=1, rowspan=1)
            self.drawRayleighDistribution(i, canv_RayleighDistribution, aH, subWaveHeightList[i - drawIndex])

            canv_FFT = plt.subplot2grid((2, 7), (1, 6), colspan=1, rowspan=1)
            self.drawFFT(i, canv_FFT)

            plt.rcParams['font.sans-serif'] = ['SimHei']
            plt.rcParams['axes.unicode_minus'] = False
            plt.subplots_adjust(left=0.03, right=0.95, top=0.95, wspace=0.212, hspace=0.45,
                                bottom=0.4)
        plt.show()

    def btn_StartProc_sy3_click(self):
        plt.close('all')
        self.ui.sy3_WaveFormOut.setText('')
        self.ui.sy3_PressFormOut.setText('')
        if self.ui.checkBox_sy3_bDrawPressForm.isChecked():
            if not self.ui.checkBox_sy3_ch_P_1.isChecked():
                self.ui.checkBox_sy3_ch_P_1.toggle()
            if not self.ui.checkBox_sy3_ch_P_2.isChecked():
                self.ui.checkBox_sy3_ch_P_2.toggle()
            if not self.ui.checkBox_sy3_ch_P_3.isChecked():
                self.ui.checkBox_sy3_ch_P_3.toggle()
            if not self.ui.checkBox_sy3_ch_P_4.isChecked():
                self.ui.checkBox_sy3_ch_P_4.toggle()
            if not self.ui.checkBox_sy3_ch_P_5.isChecked():
                self.ui.checkBox_sy3_ch_P_5.toggle()
            if not self.ui.checkBox_sy3_ch_P_6.isChecked():
                self.ui.checkBox_sy3_ch_P_6.toggle()
        self.sy3_H_rawdata = []
        self.sy3_P_rawdata = []
        for root, dirs, files in os.walk("datas\第3次海动实验\\" + self.ui.filelist_sy3.currentText() + '\\波高测量\\'):
            if files:
                for sig_file in files:
                    with open(root + sig_file, encoding='utf-8') as file:
                        sy3rawdata = [[adata.strip() for adata in data.split('    ')] for data in
                                      [data.strip() for data in file.readlines()]]
                    del (sy3rawdata[0:2])
                    del (sy3rawdata[-2:])
                    self.sy3_H_rawdata.append(sy3rawdata)

        self.sy3_H_coldatas = []
        for i_data in self.sy3_H_rawdata:
            sy3coldata = [[float(data[i]) for data in i_data] for i in range(len(i_data[0]))]
            self.sy3_H_coldatas.append(sy3coldata)
        self.sy3_timeline = list(self.floatRange(0, len(i_data), 1, 2))

        self.sy3_H_smoothdatas = []
        for sy3_t in self.sy3_H_coldatas:
            sy3_H_smooth_t = []
            for sy3_j in sy3_t:
                sy3_H_smooth_t.append(self.moving_average(sy3_j, self.ui.spinBox_smoothWindowSize_sy3.value()))
            self.sy3_H_smoothdatas.append(sy3_H_smooth_t)

        roots = []
        for root, dirs, files in os.walk(
                "datas\第3次海动实验\\" + self.ui.filelist_sy3.currentText() + '\\压力测量\\'):
            if files and (root not in roots):
                roots.append(root)
                sig_list = []
                for sig_file in files:
                    if '~' in sig_file: continue
                    sheet = xlrd.open_workbook(os.getcwd() + '\\' + root + '\\' + sig_file).sheets()[0]
                    sig_list.append((sheet.col_values(0)[1:-1], [x * 100 for x in sheet.col_values(1)[1:-1]]))
                self.sy3_P_rawdata.append(sig_list)

        self.sy3_P_smoothdatas = []
        for sy3_P_t in self.sy3_P_rawdata:
            sy3_P_smooth = []
            for sy3_P_j in sy3_P_t:
                sy3_P_smooth.append(
                    (sy3_P_j[0], self.moving_average(sy3_P_j[1], self.ui.spinBox_smoothWindowSize_sy3_P.value())))
            self.sy3_P_smoothdatas.append(sy3_P_smooth)

        drawNumber = self.ui.spinBox_number_sy3.value()
        drawIndex = self.ui.spinBox_index_sy3.value() - 1

        for i in range(drawIndex, min(drawNumber + drawIndex, len(self.sy3_H_coldatas))):
            fig = plt.figure()
            fig.canvas.set_window_title("<驻波实验> 第{}次实验".format(i + 1))
            ax_waveform = plt.subplot2grid((2, 1), (0, 0), colspan=1, rowspan=1)
            self.sy3_drawWaveform(i, ax_waveform)

            ax_Pressform = plt.subplot2grid((2, 3), (1, 0), colspan=1, rowspan=1)
            aver_P_List = self.sy3_drawPressform(i, ax_Pressform)

            ax_Featureform = plt.subplot2grid((2, 3), (1, 1), colspan=1, rowspan=1)
            self.sy3_drawPressFeatureform(i, ax_Featureform, aver_P_List)

            ax_TheoryFeatureform = plt.subplot2grid((2, 3), (1, 2), colspan=1, rowspan=1)
            self.sy3_drawPressTheoryform(i, ax_TheoryFeatureform, aver_P_List)

        plt.rcParams['font.sans-serif'] = ['SimHei']
        plt.rcParams['axes.unicode_minus'] = False
        plt.show()

    def colorSetting_click(self):
        if self.bShowColorSetting:
            self.ui.groupBox_4.hide()
            self.bShowColorSetting = False
        else:
            self.ui.groupBox_4.show()
            self.bShowColorSetting = True

    def set_WaveformColor_click(self):
        col = QColorDialog.getColor()
        if col.isValid():
            self.WaveformColor = col.name()
            self.ui.WaveformColor_btn.setStyleSheet('background: {};'.format(self.WaveformColor))

    def set_DrawPointsColor_click(self):
        col = QColorDialog.getColor()
        if col.isValid():
            self.DrawPointsColor = col.name()
            self.ui.DrawPointsColor_btn.setStyleSheet('background: {};'.format(self.DrawPointsColor))

    def set_MaxHeightPointColor_click(self):
        # self.MaxHeightPointColor = colorchooser.askcolor()[1]
        # self.ui.MaxHeightPointColor_btn.setStyleSheet('background: {};'.format(self.MaxHeightPointColor))
        col = QColorDialog.getColor()
        if col.isValid():
            self.MaxHeightPointColor = col.name()
            self.ui.MaxHeightPointColor_btn.setStyleSheet('background: {};'.format(self.MaxHeightPointColor))

    def set_MaxHeightLineColor_click(self):
        col = QColorDialog.getColor()
        if col.isValid():
            self.MaxHeightLineColor = col.name()
            self.ui.MaxHeightLineColor_btn.setStyleSheet('background: {};'.format(self.MaxHeightLineColor))

    def set_MinHeightPointColor_click(self):
        col = QColorDialog.getColor()
        if col.isValid():
            self.MinHeightPointColor = col.name()
            self.ui.MinHeightPointColor_btn.setStyleSheet('background: {};'.format(self.MinHeightPointColor))

    def set_MinHeightLineColor_click(self):
        col = QColorDialog.getColor()
        if col.isValid():
            self.MinHeightLineColor = col.name()
            self.ui.MinHeightLineColor_btn.setStyleSheet('background: {};'.format(self.MinHeightLineColor))

    def set_ZeroPointsColor_click(self):
        col = QColorDialog.getColor()
        if col.isValid():
            self.ZeroPointsColor = col.name()
            self.ui.ZeroPointsColor_btn.setStyleSheet('background: {};'.format(self.ZeroPointsColor))

    def set_datapointsColor_click(self):
        col = QColorDialog.getColor()
        if col.isValid():
            self.datapointsColor = col.name()
            self.ui.datapointsColor_btn.setStyleSheet('background: {};'.format(self.datapointsColor))

    def gif_click(self):
        col = QColorDialog.getColor()
        if col.isValid():
            self.datapointsColor = col.name()
            self.ui.datapointsColor_btn.setStyleSheet('background: {};'.format(self.datapointsColor))

    def floatRange(self, startInt, stopInt, stepInt, precision):
        f = []
        for x in range(startInt, stopInt, stepInt):
            f.append(x / (10 ** precision))
        return f

    def moving_average(self, interval, windowsize):
        window = np.ones(int(windowsize)) / float(windowsize)
        re = np.convolve(interval, window, 'same')
        return re

    def AverageSquareroot(self, subWaveHeightList):
        sumWaveH2 = 0
        for subWave in subWaveHeightList:
            sumWaveH2 += subWave ** 2
        return math.sqrt(sumWaveH2 / len(subWaveHeightList))

    def WaveFFT(self, x, t):
        fft_x = fft(x)  # fft计算
        amp_x = abs(fft_x) / len(x) * 2  # 纵坐标变换
        label_x = np.linspace(0, int(len(x) / 2) - 1, int(len(x) / 2))  # 生成频率坐标
        amp = amp_x[0:int(len(x) / 2)]  # 选取前半段计算结果即可
        # amp[0] = 0                                              # 可选择是否去除直流量信号
        fs = 1 / (t[2] - t[1])  # 计算采样频率
        fre = label_x / len(x) * fs  # 频率坐标变换
        pha = np.unwrap(np.angle(fft_x))  # 计算相位角并去除2pi跃变
        return amp, fre, pha  # 返回幅度和频率

    def SumHeigherPersent(self, subWaveHeightList, minH, maxH):
        subWaveHeightList.sort(reverse=True)
        # sum_WaveHeight = []
        # for subWave in subWaveHeightList:
        #    if(subWave >= maxH):
        #        sum_WaveHeight.append(subWave)
        #    else:
        #        return sum(sum_WaveHeight)/len(subWaveHeightList)
        num = 0
        for subWave in subWaveHeightList:
            if (subWave >= minH and subWave <= maxH):
                num += 1
        return num / len(subWaveHeightList)

    '''

    @param
    index:波高仪序号,0-5
    bRawData:True,标出滤波前的原始数据
                False:标出滤波后的平滑数据
    bDrawPoints:True:画出数据点
                False:不标出数据点
    bDrawMaxHeight:True/False 是否标出波峰
    bDrawMinHeight:True/False 是否标出波谷
    bDrawZeroPoints:True/False 是否标出零点

    '''

    def drawWaveform(self, index, canvas, bRawData, bDrawPoints, bDrawMaxHeight, bDrawMinHeight, bDrawZeroPoints):
        canvas.set_title("波形图#通道{}".format(index + 1), loc='left')
        canvas.set_xlabel("T(s)")
        canvas.set_ylabel("波高(cm)")
        ax = plt.gca()
        y_smooth = self.smooth_data[index]
        # 寻找下跨零点
        lower_crosszeropoints_X_list = []  # 下跨零点的理论X坐标
        lower_crosszeropoints_Index_list = []  # 下跨零点的数据点索引
        zeropoint_x = []
        zeropoint_y = []
        for i in range(1, len(self.timeline)):
            if (self.ui.shangkua_Radiobtn.isChecked()):
                if (y_smooth[i - 1] < 0 and y_smooth[i] > 0):
                    x_a = self.time_step * abs(y_smooth[i - 1]) / (abs(y_smooth[i - 1]) + abs(y_smooth[i]))
                    lower_crosszeropoints_Index_list.append(i)
                    zeropoint_x.append((i - 1) * self.time_step + x_a)
                    zeropoint_y.append((y_smooth[i] + y_smooth[i + 1]) * 0)
                    if len(zeropoint_x) > 1:
                        lower_crosszeropoints_X_list.append(zeropoint_x[-2:])
            else:
                if y_smooth[i - 1] > 0 and y_smooth[i] < 0:
                    x_a = self.time_step * y_smooth[i] / (y_smooth[i] + abs(y_smooth[i + 1]))
                    lower_crosszeropoints_Index_list.append(i)
                    zeropoint_x.append(i * self.time_step + x_a)
                    zeropoint_y.append((y_smooth[i] + y_smooth[i + 1]) * 0)
                    if (len(zeropoint_x) > 1):
                        lower_crosszeropoints_X_list.append(zeropoint_x[-2:])

        subWaveHeightList = []  # 每一个波周期的波高(波峰和波谷的距离H)
        subWaveHeight_Max_X_List = []
        subWaveHeight_Max_Y_List = []
        subWaveHeight_Min_X_List = []
        subWaveHeight_Min_Y_List = []

        for index_zeropoints in range(1, len(lower_crosszeropoints_Index_list)):
            subdata = list(y_smooth[
                           lower_crosszeropoints_Index_list[index_zeropoints - 1]:lower_crosszeropoints_Index_list[
                               index_zeropoints]])
            subWaveHeightList.append(abs(min(subdata)) + abs(max(subdata)))
            subWaveHeight_Max_X_List.append(
                (lower_crosszeropoints_Index_list[index_zeropoints - 1] + subdata.index(max(subdata))) * self.time_step)
            subWaveHeight_Max_Y_List.append(max(subdata))
            subWaveHeight_Min_X_List.append(
                (lower_crosszeropoints_Index_list[index_zeropoints - 1] + subdata.index(min(subdata))) * self.time_step)
            subWaveHeight_Min_Y_List.append(min(subdata))

        # print(subWaveHeightList)
        if self.ui.checkBox_bDrawWaveform.isChecked():
            if (bDrawPoints):
                if (bRawData):
                    canvas.scatter(self.timeline, self.coldata[index], s=self.ui.sy2_WaveFormPointsSize.value(),
                                   color=self.DrawPointsColor, marker='o',
                                   label="1")
                    canvas.plot(self.timeline, y_smooth, linewidth=self.ui.sy2_WaveformLineWidth.value(),
                                color=self.WaveformColor)
                else:
                    canvas.plot(self.timeline, y_smooth, marker='o', ms=self.ui.sy2_WaveFormPointsSize.value(),
                                linewidth=self.ui.sy2_WaveformLineWidth.value(), color=self.WaveformColor)
            else:
                canvas.plot(self.timeline, y_smooth,
                            linewidth=self.ui.sy2_WaveformLineWidth.value(), color=self.WaveformColor)
            if self.ui.checkBox_sy2_bWaveFormDrawLegend.isChecked():
                ax.legend(['波形'], fontsize=self.ui.sy2_waveformFontSize.value(),
                          markerscale=self.ui.sy2_MarkerSize.value())
        elif bDrawPoints:
            if (bRawData):
                canvas.scatter(self.timeline, self.coldata[index], s=self.ui.sy2_WaveFormPointsSize.value(),
                               color=self.DrawPointsColor, marker='o',
                               label="1")
                ax.legend(['数据点'], fontsize=self.ui.sy2_waveformFontSize.value(),
                          markerscale=self.ui.sy2_MarkerSize.value())
            else:
                canvas.scatter(self.timeline, y_smooth, s=self.ui.sy2_WaveFormPointsSize.value(),
                               color=self.DrawPointsColor, marker='o', label="1")
                ax.legend(['数据点'], fontsize=self.ui.sy2_waveformFontSize.value(),
                          markerscale=self.ui.sy2_MarkerSize.value())
        if bDrawMaxHeight:
            canvas.scatter(subWaveHeight_Max_X_List, subWaveHeight_Max_Y_List, s=10, color=self.MaxHeightPointColor,
                           marker='+', label="0")
            for subWaveMaxHeightIndex in range(len(subWaveHeight_Max_X_List)):
                canvas.plot([subWaveHeight_Max_X_List[subWaveMaxHeightIndex] - 1,
                             subWaveHeight_Max_X_List[subWaveMaxHeightIndex] + 1],
                            [subWaveHeight_Max_Y_List[subWaveMaxHeightIndex],
                             subWaveHeight_Max_Y_List[subWaveMaxHeightIndex]], linewidth=0.5,
                            color=self.MaxHeightLineColor)
        if bDrawMinHeight:
            canvas.scatter(subWaveHeight_Min_X_List, subWaveHeight_Min_Y_List, s=10, color=self.MinHeightPointColor,
                           marker='+', label="0")
            for subWaveMinHeightIndex in range(len(subWaveHeight_Min_X_List)):
                canvas.plot([subWaveHeight_Min_X_List[subWaveMinHeightIndex] - 1,
                             subWaveHeight_Min_X_List[subWaveMinHeightIndex] + 1],
                            [subWaveHeight_Min_Y_List[subWaveMinHeightIndex],
                             subWaveHeight_Min_Y_List[subWaveMinHeightIndex]], linewidth=0.5,
                            color=self.MinHeightLineColor)

        if bDrawZeroPoints:
            canvas.scatter(zeropoint_x, zeropoint_y, s=10, color=self.ZeroPointsColor, marker='*', label="0")

        canvas.plot([-2, 184], [0, 0], linewidth=0.5, color='blue')
        return subWaveHeightList

    def drawWaveHeightStatistics(self, index, canvas, subWaveHeightList):
        subWaveHeightList.sort(reverse=True)
        iWaveNumber = list(range(len(subWaveHeightList)))
        ax = plt.gca()
        canvas.scatter(iWaveNumber, subWaveHeightList, s=self.ui.spinBox_waveHeightPointsSize.value(), color='black',
                       marker='o', label="0")
        if self.ui.checkBox_sy2_bWaveNumberDrawLegend.isChecked():
            ax.legend(['波高'], fontsize=self.ui.sy2_wavenumberFontSize.value(),
                      markerscale=self.ui.sy2_wavenumberMarkerSize.value())
        canvas.set_title("波高统计#仪器{}".format(index + 1), loc='left')
        canvas.set_xlabel("Number of waves")
        canvas.set_ylabel("波高(cm)")
        AverageWaveHeight = sum(subWaveHeightList) / len(subWaveHeightList)
        EffectiveWaveHeight = sum(subWaveHeightList[0:math.floor(len(subWaveHeightList) / 3)]) / math.floor(
            len(subWaveHeightList) / 3)
        First10WaveHeight = sum(subWaveHeightList[0:math.floor(len(subWaveHeightList) / 10)]) / math.floor(
            len(subWaveHeightList) / 10)
        First13WaveHeight = sum(subWaveHeightList[0:math.floor(len(subWaveHeightList) * 0.13)]) / math.floor(
            len(subWaveHeightList) * 0.13)
        AverageWaveHeightSquareroot = self.AverageSquareroot(subWaveHeightList)
        if not self.ui.checkBox_drawResoult.isChecked():
            return AverageWaveHeight
        Waveout = r"平均波高:{0:.2f}cm".format(AverageWaveHeight) + '\n' + \
                  r"最大波高:{0:.2f}cm".format(max(subWaveHeightList)) + '\n' + \
                  r"有效波高:{0:.2f}cm".format(EffectiveWaveHeight) + '\n' + \
                  r"1/10波高:{0:.2f}cm".format(First10WaveHeight) + '\n' + \
                  r"13%波高:{0:.2f}cm".format(First13WaveHeight) + '\n' + \
                  r"均方根波高:{0:.2f}cm".format(AverageWaveHeightSquareroot) + '\n' + \
                  r"详情看图"

        self.ui.textEdit_WaveOut.setText(Waveout)
        canvas.text(140, 12, r"平均波高H$^-$:{0:.2f}cm".format(AverageWaveHeight))
        canvas.text(140, 10.5, r"最大波高H$_m$$_a$$_x$:{0:.2f}cm".format(max(subWaveHeightList)))
        canvas.text(140, 9, r"有效波高H$_s$:{0:.2f}cm".format(EffectiveWaveHeight))
        canvas.text(140, 7.5, r"1/10波高H$_1$$_/$$_1$$_0$:{0:.2f}cm".format(First10WaveHeight))
        canvas.text(140, 6, r"13%波高H$_1$$_3$$_\%$:{0:.2f}cm".format(First13WaveHeight))
        canvas.text(140, 4.5, r"均方根波高H$_r$$_m$$_s$:{0:.2f}cm".format(AverageWaveHeightSquareroot))

        canvas.text(140, 3,
                    r"H$_1$$_/$$_1$$_0$实验:H$_1$$_/$$_1$$_0$={0:.2f}H$^-$|理论:H$_1$$_/$$_1$$_0$=2.03H$^-$".format(
                        First10WaveHeight / AverageWaveHeight))
        canvas.text(140, 1.5, r"H$_s$实验:H$_s$={0:.2f}H$^-$|理论:H$_s$=1.60H$^-$".format(
            EffectiveWaveHeight / AverageWaveHeight))
        canvas.text(140, 0, r"H$_1$$_3$$_\%$实验:H$_1$$_3$$_\%$={0:.2f}H$^-$|理论:H$_1$$_3$$_\%$=H$^-$".format(
            First13WaveHeight / AverageWaveHeight))
        canvas.text(140, -1.5, r"H$_r$$_m$$_s$实验:H$_r$$_m$$_s$={0:.2f}H$^-$|理论:H$_r$$_m$$_s$=1.07H$^-$".format(
            AverageWaveHeightSquareroot / AverageWaveHeight))
        canvas.text(140, -3, r"H$_m$$_a$$_x$实验:H$_m$$_a$$_x$={0:.2f}H$^-$|理论:H$_m$$_a$$_x$={1:.2f}H$^-$".format(
            max(subWaveHeightList) / AverageWaveHeight, 1.07 * math.sqrt(math.log10(len(subWaveHeightList)))))
        return AverageWaveHeight

    def drawRayleighDistribution(self, index, canvas, aH, subWaveHeightList):
        canvas.set_title("瑞利分布#仪器{}".format(index + 1), loc='left')
        canvas.set_xlabel("波高(cm)")
        canvas.set_ylabel("百分比")
        ax = plt.gca()
        Rayleigh_step = self.ui.spinBox_Rayleighstep.value()
        Rayleigh_X = list(self.floatRange(0, 200, 1, 1))
        Rayleigh_Y = []
        WavePersent = []
        for i in range(0, len(Rayleigh_X), Rayleigh_step):
            Rayleigh_Y.append(
                math.pi * Rayleigh_X[i] * math.exp(-math.pi * 0.25 * (Rayleigh_X[i] / aH) ** 2) * 0.5 / aH ** 2)
            WavePersent.append(1 * self.SumHeigherPersent(subWaveHeightList, Rayleigh_X[i] - 0.2, Rayleigh_X[i] + 0.2))
        smooth_y = self.moving_average(WavePersent, self.ui.spinBox_smoothWindowSize_Rayleigh.value())
        canvas.plot(Rayleigh_X[::Rayleigh_step], Rayleigh_Y, linewidth=self.ui.sy2_RayleighLineWidth.value(),
                    color='black')
        canvas.scatter(Rayleigh_X[::Rayleigh_step], smooth_y, s=self.ui.spinBox_RayleighPointsSize.value(),
                       color=self.datapointsColor, marker='o', label="1")
        if self.ui.checkBox_sy2_bRayleighDrawLegend.isChecked():
            ax.legend(['理论曲线', '实验数据'], fontsize=self.ui.sy2_RayleighFontSize.value(),
                      markerscale=self.ui.sy2_RayleighMarkerSize.value())
        # print(Rayleigh_X[::Rayleigh_step])

    def drawFFT(self, index, canvas):
        canvas.set_title("波谱#仪器{}".format(index + 1), loc='left')
        canvas.set_xlabel("f(Hz)")
        canvas.set_ylabel("百分比")
        ax = plt.gca()
        amp, fre, pha = self.WaveFFT(self.smooth_data[index], self.timeline)  # 调用函数
        smoothed = self.moving_average(amp, self.ui.spinBox_smoothWindowSize_WaveFFT.value())
        canvas.plot(fre, smoothed, linewidth=self.ui.sy2_FFTLineWidth.value(), color='black')
        if self.ui.checkBox_sy2_bFFTDrawLegend.isChecked():
            ax.legend(['波浪谱'], fontsize=self.ui.sy2_FFTFontSize.value(),
                      markerscale=self.ui.sy2_FFTMarkerSize.value())
        Frequency_X = list(self.floatRange(0, 40, 1, 1))
        WavePu_Y = []
        for i in range(len(Frequency_X)):
            WavePu_Y.append(math.cos(Frequency_X[i] * 10) * math.sin(Frequency_X[i] * 50))
        # amp1,fre1,pha1= WaveFFT(WavePu_Y,Frequency_X)
        # canvas.plot(Frequency_X,WavePu_Y,linewidth=0.5,color='black')
        canvas.set_xlim(0, 4)
        canvas.set_ylim(0, 1)

    def sy1_drawWaveform(self, index, canvas, bUseDenoiseData, bDrawPoints):
        canvas.set_title("波形图#第{}次实验".format(index + 1), loc='left')
        canvas.set_xlabel("times(s)")
        canvas.set_ylabel("η(cm)")
        colors = ['red', 'blue', 'green', 'black', 'yellow', 'm']

        sy1_chanels = []
        if self.ui.checkBox_sy1_ch_1.isChecked():
            sy1_chanels.append(0)
        if self.ui.checkBox_sy1_ch_2.isChecked():
            sy1_chanels.append(1)
        if self.ui.checkBox_sy1_ch_3.isChecked():
            sy1_chanels.append(2)
        if self.ui.checkBox_sy1_ch_4.isChecked():
            sy1_chanels.append(3)
        if self.ui.checkBox_sy1_ch_5.isChecked():
            sy1_chanels.append(4)
        if self.ui.checkBox_sy1_ch_6.isChecked():
            sy1_chanels.append(5)
        ax = plt.gca()
        zeropoint_xs = []
        zeropoint_ys = []
        legend_line = []
        legend_list = []
        wave_heights = []
        wave_periods = []
        for i in sy1_chanels:
            y_smooth = self.sy1_smoothdatas[index][i]
            if bDrawPoints:
                if bUseDenoiseData:
                    # canvas.scatter(self.sy1_timeline, y_smooth, s=4, color=colors[i], marker='.', label="1")
                    canvas.plot(self.sy1_timeline, y_smooth, marker='o', ms=self.ui.sy1_waveformPointsScale.value(),
                                linewidth=self.ui.sy1_waveformLineWidth.value(), color=colors[i])
                    legend_line.append('波高仪#{}--波形'.format(i + 1))
                else:
                    canvas.scatter(self.sy1_timeline, self.sy1_coldatas[index][i],
                                   s=self.ui.sy1_waveformPointsScale.value(), color=colors[i], marker='o', label="1")
                    canvas.plot(self.sy1_timeline, y_smooth, linewidth=self.ui.sy1_waveformLineWidth.value(),
                                color=colors[i])
                    legend_line.append('波高仪#{}--波形'.format(i + 1))
                    legend_list.append('波高仪#{}--波高<未降噪>'.format(i + 1))
            else:
                canvas.plot(self.sy1_timeline, y_smooth, linewidth=self.ui.sy1_waveformLineWidth.value(),
                            color=colors[i])
                legend_line.append('波高仪#{}--波形'.format(i + 1))
            # canvas.plot(self.sy1_timeline, y_smooth, marker='o', ms=2, linewidth=0.5, color=colors[i])
            # 寻找下跨零点
            crosszeropoints_X_list = []  # 下跨零点的理论X坐标
            crosszeropoints_Index_list = []  # 下跨零点的数据点索引
            zeropoint_x = []
            zeropoint_y = []
            for j in range(1, len(self.sy1_timeline)):
                if self.ui.shangkua_Radiobtn_sy1.isChecked():
                    if y_smooth[j - 1] < 0 and y_smooth[j] > 0:
                        x_a = self.time_step * abs(y_smooth[j - 1]) / (abs(y_smooth[j - 1]) + y_smooth[j])
                        crosszeropoints_Index_list.append(j)
                        zeropoint_x.append((j - 1) * self.time_step + x_a)
                        zeropoint_y.append((y_smooth[j] + y_smooth[j + 1]) * 0)
                        if len(zeropoint_x) > 1:
                            crosszeropoints_X_list.append(zeropoint_x[-2:])
                else:
                    if y_smooth[j - 1] > 0 and y_smooth[j] < 0:
                        x_a = self.time_step * y_smooth[j] / (y_smooth[j] + abs(y_smooth[j + 1]))
                        crosszeropoints_Index_list.append(j)
                        zeropoint_x.append(j * self.time_step + x_a)
                        zeropoint_y.append((y_smooth[j] + y_smooth[j + 1]) * 0)
                        if len(zeropoint_x) > 1:
                            crosszeropoints_X_list.append(zeropoint_x[-2:])
            zeropoint_xs.append(zeropoint_x)
            zeropoint_ys.append(zeropoint_y)

            sy1_subWaveHeightList = []  # 每一个波周期的波高(波峰和波谷的距离H)
            for index_zeropoints in range(1, len(crosszeropoints_Index_list)):
                subdata = list(y_smooth[
                               crosszeropoints_Index_list[index_zeropoints - 1]:crosszeropoints_Index_list[
                                   index_zeropoints]])
                sy1_subWaveHeightList.append(abs(min(subdata)) + abs(max(subdata)))
            if self.ui.checkBox_sy1_bDrawZeroPoints.isChecked():
                canvas.scatter(zeropoint_x, zeropoint_y, s=self.ui.sy1_waveformZeroPointsScale.value(), color='red',
                               marker='^', label="0")
                legend_list.append('波高仪#{}--零点'.format(i + 1))
            if self.ui.checkBox_sy1_bDrawLegend.isChecked():
                ax.legend(legend_line + ['零点'], fontsize=self.ui.sy1_waveformFontSize.value(),
                          markerscale=self.ui.sy1_MarkerSize.value())

            averageWaveHeight = sum(sy1_subWaveHeightList) / len(sy1_subWaveHeightList)
            # self.ui.sy1_WaveFormOut.setText(self.ui.sy1_WaveFormOut.toPlainText()+'实验#{0}>波高传感器#{1}平均波高：{2:.2f}cm\n'.format(index+1,i+1, averageWaveHeight))
            sum_T = 0
            for zeropoint_index in range(1, len(crosszeropoints_Index_list)):
                sum_T += (crosszeropoints_Index_list[zeropoint_index] - crosszeropoints_Index_list[zeropoint_index - 1])

            self.ui.sy1_WaveFormOut.setText(
                self.ui.sy1_WaveFormOut.toPlainText() + '实验#{0}>波高传感器#{1}平均波高：{2:.2f}cm,平均波周期：{3:.4f}s\n'.format(
                    index + 1, i + 1, averageWaveHeight,
                    sum_T * self.time_step / (len(crosszeropoints_Index_list) - 1)))
            wave_heights.append(averageWaveHeight)
            wave_periods.append(sum_T * self.time_step / (len(crosszeropoints_Index_list) - 1))
        if self.ui.checkBox_sy1_bDrawAxis.isChecked():
            canvas.plot([-2, 44], [0, 0], linewidth=1, color='black')
            # print(crosszeropoints_Index_list)
        # print(zeropoint_xs)
        # 计算波速
        instrument_spacings = [0.51, 0.50, 0.61, 0.61, 0.535]
        wave_intervals = []  # 两个完整波之间的时间间隔,单位秒
        if len(zeropoint_xs) > 1:

            for zeropoint_xs_index in range(1, len(zeropoint_xs)):
                wave_interval = []
                front_wave = zeropoint_xs[zeropoint_xs_index - 1]
                back_wave = zeropoint_xs[zeropoint_xs_index]
                for zeropoint_x_front in front_wave:
                    bFindNext = False
                    for zeropoint_x_back in back_wave:
                        if bFindNext:
                            continue
                        if zeropoint_x_back > zeropoint_x_front:
                            wave_interval.append(zeropoint_x_back - zeropoint_x_front)
                            bFindNext = True
                            continue
                wave_intervals.append(wave_interval)
        self.ui.sy1_WaveFormOut.setText(
            self.ui.sy1_WaveFormOut.toPlainText() + '实验#{0}>平均波高:{1:.2f}cm,平均波周期:{2:.4f}s\n----------------------------------\n'.format(
                index + 1, sum(wave_heights) / len(wave_heights), sum(wave_periods) / len(wave_periods)))
        wave_speeds = []
        for wave_space_index in range(len(wave_intervals)):
            instrument_spacing = sum(
                instrument_spacings[sy1_chanels[wave_space_index]:sy1_chanels[wave_space_index + 1]])
            wave_spacing_time = sum(wave_intervals[wave_space_index]) / len(wave_intervals[wave_space_index])
            # print('时间:{0}s,间距:{1}m,速度:{2}m/s'.format(wave_spacing_time, instrument_spacing,instrument_spacing / wave_spacing_time))
            self.ui.sy1_WaveFormOut.setText(
                self.ui.sy1_WaveFormOut.toPlainText() + '实验#{0}>使用第#{1}波高仪与第#{2}波高仪之间dt={3:.4f}s,两仪器间隔:{4:.3f}m,\n      计算出的波速为:{5:.4f}m/s\n'.format(
                    index + 1, sy1_chanels[wave_space_index] + 1, sy1_chanels[wave_space_index + 1] + 1,
                    wave_spacing_time, instrument_spacing, instrument_spacing / wave_spacing_time))
            wave_speeds.append(instrument_spacing / wave_spacing_time)
        if len(wave_speeds) > 0:
            self.ui.sy1_WaveFormOut.setText(
                self.ui.sy1_WaveFormOut.toPlainText() + '-------------平均波速:{}m/s-----------\n'.format(
                    sum(wave_speeds) / len(wave_speeds)))

        # (2*math.pi/(sum(wave_periods) / len(wave_periods)))/((2*math.pi/(sum(wave_periods) / len(wave_periods)))**2/9.81)
        # 理论比实际速度慢?
        theorety_c = (9.81 * (sum(wave_periods) / len(wave_periods)) / 2 / math.pi) * math.tanh(
            (2 * math.pi / (sum(wave_periods) / len(wave_periods))) ** 2 * self.sy1_H / 9.81)
        theorety_L = theorety_c * (sum(wave_periods) / len(wave_periods))
        theorety_sigma = 9.81 * ((2 * math.pi / (sum(wave_periods) / len(wave_periods))) ** 2 / 9.81) * math.tanh(
            (2 * math.pi / (sum(wave_periods) / len(wave_periods))) ** 2 * self.sy1_H / 9.81)
        self.sy1_w = 2 * math.pi / (sum(wave_periods) / len(wave_periods))
        self.sy1_k = (2 * math.pi / (sum(wave_periods) / len(wave_periods))) ** 2 / 9.81
        self.ui.sy1_WaveFormOut.setText(
            self.ui.sy1_WaveFormOut.toPlainText() + '-------------Airy理论速度c = gT*tanhkh/2π = {0:.4f}m/s\n-------------Airy理论波长L = gT²*tanhkh/2π = {1:.4f}m\n-------------Airy理论波浪圆频率σ = gktanhkh = {2:.4f}\n'.format(
                theorety_c, theorety_L, theorety_sigma))
        self.ui.sy1_WaveFormOut.setText(
            self.ui.sy1_WaveFormOut.toPlainText() + '实验#{0}👆>--------------------------end-------------------\n'.format(
                index + 1))
        # canvas.set_xlim(15,20)
        canvas.set_ylim(-max(sy1_subWaveHeightList) * 0.5 - self.ui.sy1_waveformYScale.value(),
                        max(sy1_subWaveHeightList) * 0.5 + self.ui.sy1_waveformYScale.value())

    def nStokesWave(self, n, A, x, t, w, k):
        ad = []
        for index in range(1, n):
            ad.append((A ** index * k ** (index - 1) * math.cos(index * (k * x - w * t)) / index))
        return sum(ad)

    def sy1_ContrastWaveForm(self, index, canvas):
        canvas.set_title("波形图#第{}次实验".format(index + 1), loc='left')
        canvas.set_xlabel("times(s)")
        canvas.set_ylabel("η(cm)")

        ch = 0
        if self.ui.sy1_contract_ch_1.isChecked():
            ch = 0
        if self.ui.sy1_contract_ch_2.isChecked():
            ch = 1
        if self.ui.sy1_contract_ch_3.isChecked():
            ch = 2
        if self.ui.sy1_contract_ch_4.isChecked():
            ch = 3
        if self.ui.sy1_contract_ch_5.isChecked():
            ch = 4
        if self.ui.sy1_contract_ch_6.isChecked():
            ch = 5

        ax = plt.gca()
        zeropoint_xs = []
        zeropoint_ys = []
        legend_line = []
        legend_list = []
        wave_heights = []
        wave_periods = []

        y_smooth = self.sy1_smoothdatas[index][ch]
        canvas.scatter(self.sy1_timeline, y_smooth, s=self.ui.sy1_constrastPointsScale.value(), color='black',
                       marker='o', label="1")
        # canvas.plot(self.sy1_timeline, y_smooth, marker='o', ms=self.ui.sy1_waveformPointsScale.value(),
        #            linewidth=self.ui.sy1_waveformLineWidth.value(), color='black')

        # 寻找零点
        crosszeropoints_X_list = []  # 下跨零点的理论X坐标
        crosszeropoints_Index_list = []  # 下跨零点的数据点索引
        zeropoint_x = []
        zeropoint_y = []
        for j in range(1, len(self.sy1_timeline)):
            if self.ui.shangkua_Radiobtn_sy1.isChecked():
                if y_smooth[j - 1] < 0 < y_smooth[j]:
                    x_a = self.time_step * abs(y_smooth[j - 1]) / (abs(y_smooth[j - 1]) + y_smooth[j])
                    crosszeropoints_Index_list.append(j)
                    zeropoint_x.append((j - 1) * self.time_step + x_a)
                    zeropoint_y.append((y_smooth[j] + y_smooth[j + 1]) * 0)
                    if len(zeropoint_x) > 1:
                        crosszeropoints_X_list.append(zeropoint_x[-2:])
            else:
                if y_smooth[j - 1] > 0 > y_smooth[j]:
                    x_a = self.time_step * y_smooth[j] / (y_smooth[j] + abs(y_smooth[j + 1]))
                    crosszeropoints_Index_list.append(j)
                    zeropoint_x.append(j * self.time_step + x_a)
                    zeropoint_y.append((y_smooth[j] + y_smooth[j + 1]) * 0)
                    if len(zeropoint_x) > 1:
                        crosszeropoints_X_list.append(zeropoint_x[-2:])

        H = 2
        if index == 0:
            H = self.ui.sy1_contrastform_H_1.value()
        elif index == 1:
            H = self.ui.sy1_contrastform_H_2.value()
        elif index == 2:
            H = self.ui.sy1_contrastform_H_3.value()
        A = H / 2
        AiryWave = lambda x, t, w, k: A * math.cos(w * t - k * x)
        StokesWave = lambda x, t, w, k: 1 / 4 * A ** 2 * k * math.cos(2 * (k * x - w * t)) + A * math.cos(
            w * t - k * x)
        Airy_y = []
        Stokes_y = []
        q = 2 * math.pi / self.sy1_k
        T = 2 * math.pi / self.sy1_w
        if self.ui.shangkua_Radiobtn_sy1.isChecked():
            shift = (self.sy1_w * (
                    zeropoint_x[2] % T) - math.pi * 0.5) / self.sy1_k + 0.5 * q + self.ui.sy1_contrastform_Shift.value()
        else:
            shift = (self.sy1_w * (
                    zeropoint_x[2] % T) + math.pi * 0.5) / self.sy1_k + 0.5 * q + self.ui.sy1_contrastform_Shift.value()
        # print((self.sy1_w*(zeropoint_x[12]% T)-math.pi*0.5)/self.sy1_k + 0.25*T, 0.6*T)
        for tl in self.sy1_timeline:
            Airy_y.append(AiryWave(shift, tl, self.sy1_w, self.sy1_k))
            Stokes_y.append(StokesWave(shift, tl, self.sy1_w, self.sy1_k))
        print(index)
        canvas.plot(self.sy1_timeline, Airy_y, linewidth=self.ui.sy1_constrastLineWidth.value(),
                    color='black')
        canvas.plot(self.sy1_timeline, Stokes_y, linewidth=self.ui.sy1_constrastLineWidth.value(),
                    color='red')
        # canvas.plot([-2, 44], [0, 0], linewidth=1, color='black')
        ax.legend(['Airy', '2rd Stocks', 'Experiments'], fontsize=self.ui.sy1_waveformFontSize.value(),
                  markerscale=self.ui.sy1_MarkerSize.value())
        canvas.set_xlim(min(self.ui.sy1_contrastform_ShowRange.value(), self.ui.sy1_contrastform_ShowRange_2.value()),
                        max(self.ui.sy1_contrastform_ShowRange.value(), self.ui.sy1_contrastform_ShowRange_2.value()))

    def sy3_drawWaveform(self, index, canvas):
        canvas.set_title("波形图#第{}次实验".format(index + 1), loc='left')
        canvas.set_xlabel("times(s)")
        canvas.set_ylabel("η(cm)")
        colors = ['red', 'blue', 'green', 'black', 'yellow', 'm']

        sy3_chanels = []
        if self.ui.checkBox_sy3_ch_1.isChecked():
            sy3_chanels.append(0)
        if self.ui.checkBox_sy3_ch_2.isChecked():
            sy3_chanels.append(1)
        if self.ui.checkBox_sy3_ch_3.isChecked():
            sy3_chanels.append(2)
        if self.ui.checkBox_sy3_ch_4.isChecked():
            sy3_chanels.append(3)
        if self.ui.checkBox_sy3_ch_5.isChecked():
            sy3_chanels.append(4)
        if self.ui.checkBox_sy3_ch_6.isChecked():
            sy3_chanels.append(5)

        ax = plt.gca()
        zeropoint_xs = []
        zeropoint_ys = []
        legend_line = []
        legend_list = []
        wave_heights = []
        wave_periods = []
        for i in sy3_chanels:
            y_smooth = self.sy3_H_smoothdatas[index][i]
            if self.ui.checkBox_sy3_bDrawPoints.isChecked():
                if self.ui.checkBox_sy3_bSmoothData.isChecked():
                    canvas.plot(self.sy3_timeline, y_smooth, marker='o', ms=self.ui.sy3_waveformPointsScale.value(),
                                linewidth=self.ui.sy3_waveformLineWidth.value(), color=colors[i])
                else:
                    canvas.scatter(self.sy3_timeline, self.sy3_H_coldatas[index][i],
                                   s=self.ui.sy3_waveformPointsScale.value(), color=colors[i], marker='o', label="1")
                    canvas.plot(self.sy3_timeline, y_smooth, linewidth=self.ui.sy3_waveformLineWidth.value(),
                                color=colors[i])
            else:
                canvas.plot(self.sy3_timeline, y_smooth, linewidth=self.ui.sy3_waveformLineWidth.value(),
                            color=colors[i])

            # 寻找零点
            crosszeropoints_X_list = []  # 下跨零点的理论X坐标
            crosszeropoints_Index_list = []  # 下跨零点的数据点索引
            zeropoint_x = []
            zeropoint_y = []
            for j in range(1, len(self.sy3_timeline)):
                if self.ui.shangkua_Radiobtn_sy3.isChecked():
                    if y_smooth[j - 1] < 0 < y_smooth[j]:
                        x_a = self.time_step * abs(y_smooth[j - 1]) / (abs(y_smooth[j - 1]) + y_smooth[j])
                        crosszeropoints_Index_list.append(j)
                        zeropoint_x.append((j - 1) * self.time_step + x_a)
                        zeropoint_y.append((y_smooth[j] + y_smooth[j + 1]) * 0)
                        if len(zeropoint_x) > 1:
                            crosszeropoints_X_list.append(zeropoint_x[-2:])
                else:
                    if y_smooth[j - 1] > 0 > y_smooth[j]:
                        x_a = self.time_step * y_smooth[j] / (y_smooth[j] + abs(y_smooth[j + 1]))
                        crosszeropoints_Index_list.append(j)
                        zeropoint_x.append(j * self.time_step + x_a)
                        zeropoint_y.append((y_smooth[j] + y_smooth[j + 1]) * 0)
                        if len(zeropoint_x) > 1:
                            crosszeropoints_X_list.append(zeropoint_x[-2:])
            zeropoint_xs.append(zeropoint_x)
            zeropoint_ys.append(zeropoint_y)

            sy3_subWaveHeightList = []  # 每一个波周期的波高(波峰和波谷的距离H)
            for index_zeropoints in range(1, len(crosszeropoints_Index_list)):
                subdata = list(y_smooth[
                               crosszeropoints_Index_list[index_zeropoints - 1]:crosszeropoints_Index_list[
                                   index_zeropoints]])
                sy3_subWaveHeightList.append(abs(min(subdata)) + abs(max(subdata)))
            if self.ui.checkBox_sy3_bDrawZeroPoints.isChecked():
                canvas.scatter(zeropoint_x, zeropoint_y, s=self.ui.sy3_waveformZeroPointsScale.value(), color='red',
                               marker='^', label="0")
            averageWaveHeight = sum(sy3_subWaveHeightList) / len(sy3_subWaveHeightList)
            # self.ui.sy1_WaveFormOut.setText(self.ui.sy1_WaveFormOut.toPlainText()+'实验#{0}>波高传感器#{1}平均波高：{2:.2f}cm\n'.format(index+1,i+1, averageWaveHeight))
            sum_T = 0
            for zeropoint_index in range(1, len(crosszeropoints_Index_list)):
                sum_T += (crosszeropoints_Index_list[zeropoint_index] - crosszeropoints_Index_list[zeropoint_index - 1])

            self.ui.sy3_WaveFormOut.setText(
                self.ui.sy3_WaveFormOut.toPlainText() + '实验#{0}>波高传感器#{1}平均波高：{2:.2f}cm,平均波周期：{3:.4f}s\n'.format(
                    index + 1, i + 1, averageWaveHeight,
                    sum_T * self.time_step / (len(crosszeropoints_Index_list) - 1)))

    def sy3_drawPressform(self, index, canvas):
        canvas.set_title("波压力图#第{}次实验".format(index + 1), loc='left')
        canvas.set_xlabel("times(s)")
        canvas.set_ylabel("P(Kpa)")
        colors = ['red', 'blue', 'green', 'black', 'yellow', 'm']

        sy3_chanels = []
        if self.ui.checkBox_sy3_ch_P_1.isChecked():
            sy3_chanels.append(0)
        if self.ui.checkBox_sy3_ch_P_2.isChecked():
            sy3_chanels.append(1)
        if self.ui.checkBox_sy3_ch_P_3.isChecked():
            sy3_chanels.append(2)
        if self.ui.checkBox_sy3_ch_P_4.isChecked():
            sy3_chanels.append(3)
        if self.ui.checkBox_sy3_ch_P_5.isChecked():
            sy3_chanels.append(4)
        if self.ui.checkBox_sy3_ch_P_6.isChecked():
            sy3_chanels.append(5)

        averagePressList = []
        for i in sy3_chanels:
            x_axis = self.sy3_P_smoothdatas[index][i][0]
            y_smooth = self.sy3_P_smoothdatas[index][i][1]
            canvas.plot(x_axis, y_smooth, marker='o', ms=self.ui.sy3_waveformPointsScale_P.value(),
                        linewidth=self.ui.sy3_waveformLineWidth_P.value(), color=colors[i])
            y_shift = sum(y_smooth) / len(y_smooth) + self.ui.sy3_PressformShift.value()
            if self.ui.checkBox_sy3_bDrawDivideLine.isChecked():
                canvas.plot([x_axis[0], x_axis[-1]], [y_shift, y_shift], linewidth=1, color='black')

            crosszeropoints_Index_list = []  # 下跨零点的数据点索引
            zeropoint_x = []
            zeropoint_y = []
            for j in range(1, len(x_axis)):
                if y_smooth[j - 1] < y_shift < y_smooth[j]:
                    crosszeropoints_Index_list.append(j)
                    zeropoint_x.append(x_axis[j])
                    zeropoint_y.append(y_shift)

            sy3_subPressList = []  # 每一个周期的最大波压力(Kpa)
            maxpoint_x = []
            maxpoint_y = []
            for index_zeropoints in range(1, len(crosszeropoints_Index_list)):
                subdata = list(y_smooth[
                               crosszeropoints_Index_list[index_zeropoints - 1]:crosszeropoints_Index_list[
                                   index_zeropoints]])
                n_index = crosszeropoints_Index_list[index_zeropoints - 1] + subdata.index(max(subdata))
                maxpoint_x.append(x_axis[n_index])
                maxpoint_y.append(y_smooth[n_index])
                sy3_subPressList.append(max(subdata))
            averagePress = 0
            if len(sy3_subPressList) > 0:
                averagePress = sum(sy3_subPressList) / len(sy3_subPressList)
            averagePressList.append(averagePress)

            if self.ui.checkBox_sy3_bDrawDivideLine_2.isChecked():
                canvas.scatter(zeropoint_x, zeropoint_y, s=self.ui.sy3_waveformZeroPointsScale.value(), color='red',
                               marker='^', label="0")
            if self.ui.checkBox_sy3_bDrawDivideLine_3.isChecked():
                canvas.scatter(maxpoint_x, maxpoint_y, s=self.ui.sy3_waveformZeroPointsScale.value(), color='red',
                               marker='^', label="0")

            self.ui.sy3_PressFormOut.setText(
                self.ui.sy3_PressFormOut.toPlainText() + '实验#{0}>压力传感器#{1}平均压力：{2:.4f}Kpa\n'.format(
                    index + 1, i + 1, averagePress))
        self.ui.sy3_PressFormOut.setText(
            self.ui.sy3_PressFormOut.toPlainText() + '实验#{0}>最大压力在#{1}号传感器处,为:{2:.4f}Kpa\n'.format(
                index + 1, averagePressList.index(max(averagePressList)) + 1, max(averagePressList)))
        self.ui.sy3_PressFormOut.setText(
            self.ui.sy3_PressFormOut.toPlainText() + '实验#{0}>↑-----------------------------end-------------↑\n'.format(
                index + 1))
        return averagePressList

    def sy3_drawPressFeatureform(self, index, canvas, aver_P_List):
        canvas.set_title("波压力特征图#第{}次实验".format(index + 1), loc='left')
        canvas.set_xlabel("P(Kpa)")
        canvas.set_ylabel("H(cm)")
        colors = ['red', 'blue', 'green', 'black', 'yellow', 'm']

        # 绘制直立堤
        X_Shift = self.ui.sy3_PressForm_X_Shift.value()
        LineWidth = self.ui.sy3_PressFormLineWidth.value()

        canvas.plot([0, 80], [0, 0], linewidth=LineWidth, color='black')
        canvas.plot([50 + X_Shift, 50 + X_Shift, 55 + X_Shift, 60 + X_Shift, 80 + X_Shift, 80 + X_Shift],
                    [0, 60, 60, 50, 50, 0], linewidth=LineWidth, color='black')
        # 绘制水面
        if self.ui.checkBox_sy3_bDrawWater.isChecked():
            canvas.plot([0, 50 + X_Shift], [30, 30], linewidth=LineWidth, color='blue')

        # 绘制波压力特征图
        PressFeature_X = [(50 + X_Shift) - x * self.ui.sy3_PressForm_Scale.value() for x in aver_P_List]
        PressFeature_Y = [16, 21, 26, 31, 36, 41]

        canvas.plot(PressFeature_X, PressFeature_Y, linewidth=LineWidth, color='blue')

    def sy3_drawPressTheoryform(self, index, canvas, aver_P_List):
        canvas.set_title("波压力特征图#第{}次实验".format(index + 1), loc='left')
        canvas.set_xlabel("P(kPa)")
        canvas.set_ylabel("z(cm)")
        colors = ['red', 'blue', 'green', 'black', 'yellow', 'm']

        PressTheory_X = aver_P_List[0:3]
        PressTheory_Y = [x-30 for x in [16, 21, 26, 31, 36, 41][0:3]]

        canvas.scatter(PressTheory_X, PressTheory_Y, s=self.ui.sy3_waveformZeroPointsScale.value(), color='red',
                       marker='^', label="0")

if __name__ == "__main__":
    QCoreApplication.setAttribute(Qt.AA_EnableHighDpiScaling)
    app = QApplication(sys.argv)
    app.setWindowIcon(QIcon("res\\coastal.ico"))

    gui = CoastalDataProc()
    gui.ui.show()

    app.lastWindowClosed.connect(gui.close_Event)

    gui.ui.StartProc_sy1.clicked.connect(gui.btn_StartProc_sy1_click)
    gui.ui.StartProc_sy2.clicked.connect(gui.btn_StartProc_sy2_click)
    gui.ui.StartProc_sy3.clicked.connect(gui.btn_StartProc_sy3_click)

    gui.ui.WaveformColor_btn.clicked.connect(gui.set_WaveformColor_click)
    gui.ui.DrawPointsColor_btn.clicked.connect(gui.set_DrawPointsColor_click)
    gui.ui.MaxHeightPointColor_btn.clicked.connect(gui.set_MaxHeightPointColor_click)
    gui.ui.MaxHeightLineColor_btn.clicked.connect(gui.set_MaxHeightLineColor_click)
    gui.ui.MinHeightPointColor_btn.clicked.connect(gui.set_MinHeightPointColor_click)
    gui.ui.MinHeightLineColor_btn.clicked.connect(gui.set_MinHeightLineColor_click)
    gui.ui.ZeroPointsColor_btn.clicked.connect(gui.set_ZeroPointsColor_click)
    gui.ui.datapointsColor_btn.clicked.connect(gui.set_datapointsColor_click)
    gui.ui.show_colorsetting.clicked.connect(gui.colorSetting_click)
    sys.exit(app.exec_())
    # app.exec_()
