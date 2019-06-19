#!/usr/bin/env python3
from tkinter import *
from tkinter import ttk
from nei import *
import os
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import threading


def read_default():
    with open(Path('MU\materials\default.txt'), 'r') as file:
        content = file.read().upper()
    materials = re.split('\s*\n\s*', content)
    return materials


def read_materials_file(filename,verbose=False):
    with open(Path('MU\\materials\\' + filename + '.txt'), 'r') as file:
        content = file.read().upper()
    materials = re.split('\s*\n\s*', content)
    return materials


def input_materials():
    materials = input('\nPlease input the names of the materials to investigate.\n'
                      'For example: K2SeO4 Se-Meth, Water\n'
                      'Or press Enter to skip\n')
    materials = re.findall(r"[\w'-]+", materials)
    return materials


def write_materials_file(materials, filename):
    filename = filename.lower()
    with open(Path('MU\\materials\\' + filename + '.txt'), 'w') as file:
        if type(materials) == 'str':
            file.write(materials)
        else:
            file.write('\n'.join(materials))
    return


class gui_get_materials:
    """docstring"""

    def __init__(self, save_path=''):
        #################### create some ttk styles     ######################
        # s = ttk.Style()
        # s.configure("Small.TNotebook.Tab", padding=[25,5])
        # # s.configure("BW.TNotebook.Tab", padding=[157,5])
        # # s.configure('Kim.TNotebook.Tab', padding=[25,5])
        #################### Start window   ##################################
        self.save_path = save_path
        # self.window = Tk()
        self.window = Toplevel()
        self.window.minsize(width=400, height=300)
        self.window.title("Input the materials")
        master = Frame(self.window, width=300)
        # master = ttk.Notebook(window)
        # master.pack(expand=1, fill="both")
        master.pack()
        ######################################################################
        #################### frame input field ###################
        # Main frame
        Label(master, text='CHOOSE FROM EXISTING MATERIALS FILE').pack()
        frame0 = Frame(master, bd=4)
        frame0.pack(pady=5)
        btSelectFile = Button(frame0, text='Browse', command=self.select)
        btSelectFile.grid(row=0, columnspan=2)

        Label(master, text='MATERIALS.\n(Hit ENTER on the keyboard after every material)').pack()
        frame1 = Frame(master, bd=4)
        frame1.pack(pady=5)
        sbText = Scrollbar(frame1)
        self.textMaterials = Text(frame1, height=9, width=30, yscrollcommand=sbText.set,
                                  relief=SUNKEN)
        self.textMaterials.grid(row=0, column=0)
        sbText.config(command=self.textMaterials.yview)
        sbText.grid(row=0, column=1, sticky=NS)
        ######################### Frame2 for materials file name ##########
        Label(master, text='MATERIALS ALIAS').pack(pady=5)
        frame2 = Frame(master, bd=4)
        frame2.pack(pady=5)
        self.filename = Entry(frame2, width=30)
        self.filename.pack()
        ######################### FRAME3 for others #######################
        # Main frame
        frame3 = Frame(master, bd=4)
        frame3.pack(pady=10)
        # frame3.grid(row=2, column=0, sticky=NS)
        # TODO: frame for REMEMBER LAST SETTING
        btConfirm = Button(frame3, text='CONFIRM', command=self.confirm)
        Button.config(btConfirm, font=('Helvetica', '11'))
        btConfirm.grid(row=0, columnspan=3)
        master.mainloop()

    def select(self):
        filename = filedialog.askopenfilename(title='Please select file')
        print('Selected materials file:',filename)
        # import os
        filename = os.path.basename(filename).split('.')[0]  # get the filename from the full path
        materials = read_materials_file(filename)
        write_materials_file(materials, 'last')
        self.window.destroy()

    def confirm(self):
        materials = self.textMaterials.get('1.0',
                                           'end-1c')  # https://stackoverflow.com/questions/14824163/how-to-get-the-input-from-the-tkinter-text-box-widget
        materials = re.split('\s*\n\s*', materials)
        filename = self.filename.get()
        write_materials_file(materials, filename)
        write_materials_file(materials, filename='last')
        self.window.destroy()


class gui_get_arrangement():
    """docstring"""

    def __init__(self, root, save_path=''):
        #################### create some ttk styles     ######################
        # s = ttk.Style()
        # s.configure("Small.TNotebook.Tab", padding=[25,5])
        # # s.configure("BW.TNotebook.Tab", padding=[157,5])
        # # s.configure('Kim.TNotebook.Tab', padding=[25,5])

        #################### Start window   ##################################

        self.save_path = save_path
        self.window = Toplevel()
        self.window.minsize(width=400, height=450)
        self.window.title("SYSTEM & DETECTOR ARRANGEMENTS")
        master = Frame(self.window, width=300)
        # master = ttk.Notebook(window)
        # master.pack(expand=1, fill="both")
        master.pack()
        ######################################################################

        #################### frame1 for system arrangement ###################
        # Main frame
        # frame1 = LabelFrame(master,width=500, text = 'SYSTEM',bd=4, relief=GROOVE,labelanchor=N)
        Label(master, text='SYSTEM').pack(pady=5)

        frame1 = Frame(master, width=700)
        # frame1.grid(row=0, column=0, sticky=NS)
        frame1.pack(pady=5)
        # Label(frame1, text='').grid()

        # frame for the ARRANGEMENT NAME.
        frame1_1 = Frame(frame1)
        frame1_1.grid(sticky=W)
        Label(frame1_1, text='NAME', width=20).grid(row=0, sticky=W)

        self.aName = StringVar()
        # self.aName.set('Default')
        self.entryName = Entry(frame1_1, textvariable=self.aName)
        self.entryName.grid(row=0, column=1)
        self.aName.set('Default')

        # frame for DIFFRACTION PLANE
        frame1_2 = Frame(frame1)
        frame1_2.grid(sticky=W)
        Label(frame1_2, text='DIFFRACTION PLANE', width=20).grid(row=0, sticky=W)

        self.diffPlane = StringVar()
        self.rbtVerti = Radiobutton(frame1_2, text='Vertical', variable=self.diffPlane, value='Vertical')
        self.rbtVerti.grid(row=0, column=1)

        self.rbtHoriz = Radiobutton(frame1_2, text='Horizontal', variable=self.diffPlane, value='Horizontal')
        self.rbtHoriz.grid(row=0, column=2)

        self.diffPlane.set('Vertical')

        # frame for ASYMMETRY ANGLE
        frame1_3 = Frame(frame1)
        frame1_3.grid(sticky=W)
        Label(frame1_3, text='ASYMMETRY ANGLE', width=20).grid(row=0, sticky=W)

        self.chi = DoubleVar()
        self.entryChi = Entry(frame1_3, width=12, textvariable=self.chi).grid(row=0, column=1)
        self.chi.set(0.0)

        Label(frame1_3, text='Degree').grid(row=0, column=2)

        # frame for H,K,L
        frame1_4 = Frame(frame1)
        frame1_4.grid(sticky=W)
        Label(frame1_4, text='H,K,L', width=20).grid(row=0, sticky=W)

        self.h = IntVar();
        self.k = IntVar();
        self.l = IntVar()
        Label(frame1_4, width=3, text='H').grid(row=0, column=1)
        self.entryH = Entry(frame1_4, width=3, textvariable=self.h).grid(row=0, column=2)
        self.h.set(1)

        Label(frame1_4, width=3, text='K').grid(row=0, column=3)
        self.entryK = Entry(frame1_4, width=3, textvariable=self.k).grid(row=0, column=4)
        self.k.set(1)

        Label(frame1_4, width=3, text='L').grid(row=0, column=5)
        self.entryL = Entry(frame1_4, width=3, textvariable=self.l).grid(row=0, column=6)
        self.l.set(1)

        # frame for K-EDGE
        frame1_5 = Frame(frame1)
        frame1_5.grid(sticky=W)
        Label(frame1_5, text='K-EDGE', width=20).grid(row=0, sticky=W)

        self.kEdge = DoubleVar()
        self.entryEdge = Entry(frame1_5, width=12, textvariable=self.kEdge).grid(row=0, column=1)
        self.kEdge.set(0.0)

        Label(frame1_5, text='keV').grid(row=0, column=2)

        # frame for ENERGY RANGE
        frame1_6 = Frame(frame1)
        frame1_6.grid(sticky=W)
        Label(frame1_6, text='ENERGY RANGE', width=20).grid(row=0, sticky=W)

        Label(frame1_6, text='LOW', width=5).grid(row=0, column=1)
        self.lowE = DoubleVar()
        self.entryLowE = Entry(frame1_6, width=5, textvariable=self.lowE).grid(row=0, column=2)
        self.lowE.set(0.0)

        Label(frame1_6, text='HIGH', width=5).grid(row=0, column=3)
        self.highE = DoubleVar()
        self.entryHighE = Entry(frame1_6, width=5, textvariable=self.highE).grid(row=0, column=4)
        self.highE.set(0.0)

        Label(frame1_6, text='keV').grid(row=0, column=5)

        # frame for DISTANCE between FOCUS and DETECTOR
        frame1_7 = Frame(frame1)
        frame1_7.grid(sticky=W)
        Label(frame1_7, text='FOCUS to DETECTOR', width=20).grid(row=0, sticky=W)

        self.f2d = DoubleVar()
        self.entryf2d = Entry(frame1_7, width=12, textvariable=self.f2d).grid(row=0, column=1)
        self.f2d.set(0.0)

        Label(frame1_7, text='mm').grid(row=0, column=2)

        separator = Frame(master, height=2, bd=1, relief=SUNKEN)
        separator.pack(fill=X, padx=5, pady=5)

        #################### frame2 for detector arrangement ###################
        # Main frame
        # frame2 = LabelFrame(master, width=50, text='DETECTOR', bd=4, relief=GROOVE, labelanchor=N)
        Label(master, text='DETECTOR').pack()

        frame2 = Frame(master, width=800)
        # frame2.grid(row=1, column=0, sticky=NS)
        frame2.pack(pady=10)
        # Label(frame2, text='').grid() #Make a blank line

        # frame for DETECTOR TYPE
        frame2_1 = Frame(frame2)
        frame2_1.grid(sticky=W)
        Label(frame2_1, text='TYPE', width=20).grid(row=0, sticky=W)

        self.aType = StringVar()
        self.entryType = Entry(frame2_1, textvariable=self.aType)
        self.entryType.grid(row=0, column=1)
        self.aType.set('Default')

        # frame for PIXEL SIZE
        # frame for DISTANCE between FOCUS and DETECTOR
        frame2_2 = Frame(frame2)
        frame2_2.grid(sticky=W)
        Label(frame2_2, text='PIXEL SIZE', width=20).grid(row=0, sticky=W)

        self.pixelSize = DoubleVar()
        self.entryPixelSize = Entry(frame2_2, width=12, textvariable=self.pixelSize).grid(row=0, column=1)
        self.kEdge.set(0.0)

        Label(frame2_2, text='um').grid(row=0, column=2)

        # frame for DETECTOR THRESHOLD
        frame2_3 = Frame(frame2)
        frame2_3.grid(sticky=W)
        Label(frame2_3, text='THRESHOLD (%)', width=20).grid(row=0, sticky=W)

        self.thres = DoubleVar()
        self.entryThres = Entry(frame2_3, width=12, textvariable=self.thres).grid(row=0, column=1)
        self.thres.set(50.0)

        separator = Frame(master, height=2, bd=1, relief=SUNKEN)
        separator.pack(fill=X, padx=5, pady=5)

        ######################### get the values  ####################

        # self.diffaction_plane = self.diffPlane.get()
        # self.type = self.aName.get()
        # self.chi_degrees = self.chi.get()
        # self.hkl = [self.h.get(), self.k.get(), self.l.get()]
        self.energy = self.kEdge.get()
        # self.energy_range = [self.lowE.get(),self.highE.get()]
        # self.dist_fd = self.f2d.get()

        ######################### FRAME3 for others #######################
        # Main frame
        frame3 = Frame(master, bd=4)
        frame3.pack(pady=10)
        # frame3.grid(row=2, column=0, sticky=NS)
        # TODO: frame for REMEMBER LAST SETTING

        # TODO: frame for the CONFIRM button
        btConfirm = Button(frame3, text='CONFIRM', command=self.confirm)
        Button.config(btConfirm, font=('Helvetica', '11'))
        btConfirm.grid(row=0, columnspan=3)
        master.mainloop()

    def confirm(self):
        #     self.diffaction_plane = self.diffPlane
        #     self.type = self.aName
        #     self.chi_degrees = self.chi
        #     self.hkl = [self.h, self.k, self.l]
        #     self.energy = self.kEdge
        #     self.energy_range = [self.lowE,self.highE]
        #     self.dist_fd = self.f2d
        #     # self.detector = self.detector(data)
        aDict = {}

        aDict['type'] = self.aName.get()
        aDict['diffraction_plane'] = self.diffPlane.get()
        aDict['chi_degrees'] = self.chi.get()
        aDict['h'] = self.h.get()
        aDict['k'] = self.k.get()
        aDict['l'] = self.l.get()
        aDict['energy'] = self.kEdge.get()
        aDict['energy_range_low'] = self.lowE.get()
        aDict['energy_range_high'] = self.highE.get()
        aDict['dist_fd'] = self.f2d.get()
        aDict['det_type'] = self.aType.get()
        aDict['det_pixel'] = self.pixelSize.get() / 1000
        aDict['det_pct_max'] = self.thres.get()
        # todo: Be careful. Some values for arrangement are set to zero, because it is not in the gui
        aDict['det_flip'] = 0
        aDict['det_phperapu'] = 0
        aDict['det_disp_x_demag'] = 0

        arrange_df = pd.DataFrame.from_dict(aDict, orient='index')
        arrange_df.to_csv(Path(self.save_path) / 'arrangement.dat', header=False)
        # print(arrange_df)
        self.window.destroy()

########################  Main GUI  ###########################
class NearEdgeImaging:
    """docstring for NearEdgeImaging"""

    def __init__(self):
        #################### create some ttk styles     ######################
        # s = ttk.Style()
        # s.configure("Small.TNotebook.Tab", padding=[25,5])
        # # s.configure("BW.TNotebook.Tab", padding=[157,5])
        # # s.configure('Kim.TNotebook.Tab', padding=[25,5])


        #################### Start window   ##################################
        self.window = Tk()
        self.window.minsize(width=800, height=500)
        self.window.title("Near Edge Imaging")

        master = ttk.Notebook(self.window)
        # master.pack(expand=1, fill="both")
        master.grid()
        tabKES = ttk.Frame(master)
        #################  Frame 0  ### Jun 18, 2019 Todo ######################
        frame0 = LabelFrame(tabKES,bd=4,relief=RIDGE,
                            text="STEP 0: SYSTEM ARRANGEMENT and MATERIALS",
                            font=('Helvetica', '11'), labelanchor=N)
        frame0.pack(pady=12)

        # This block is used to create some grid spacing (Begin)
        rows = 0
        while rows < 50:
            frame0.rowconfigure(rows, weight=1)
            frame0.columnconfigure(rows, weight=1)
            rows += 1
        # This block is used to create some grid spacing (End)

        btArrangement = Button(frame0, text="System Arrangement", command=self.call_arrangement_gui)
        btArrangement.grid(row=0,column=0)

        # If the materials button is not clicked, then the program will use default file.

        self.select_materials = BooleanVar()
        self.select_materials.set(False) # Material file is not specified
        btMaterials = Button(frame0, text='Materials', command=self.call_materials_gui)
        btMaterials.grid(row=0, column=6,columnspan=2)

        #
        # btSelect = Button(frame0,text='Select', command=self.select)
        # btSelect.grid(row=0, column=2)
        #################   Frame 1   ##########################################

        frame1 = LabelFrame(tabKES,bd=4,relief=RIDGE,
                            text = "STEP 1: SETUP PARAMETERS\n(Use default if you are not sure what they are)",
                            font = ('Helvetica', '11'),labelanchor=N)
        frame1.pack(pady=12)
        # # materials
        # ### checkbutton for 'set as default'
        # self.set_default = BooleanVar()
        # cbtDefault = Checkbutton(frame1, text='Set as default',variable=self.set_default)
        # self.set_default.set(False)
        #
        # ### list view for 'choose materials set'
        # self.cboxMatSets = ttk.Combobox()
        # self.cboxMatSets.current(0)#show the first element in the combobox
        # ### button to input 'new materials'

        # use file
        self.use_file = BooleanVar()
        cbtUseFile = Checkbutton(frame1, text="Use File", variable=self.use_file)
        self.use_file.set(True)

        # lowpass filter
        self.lowpass = BooleanVar()
        cbtLowpass = Checkbutton(frame1, text="Lowpass Filter", variable=self.lowpass)

        # multislice
        self.multislice = BooleanVar()
        cbtMultislice = Checkbutton(frame1, text='Multiple Slices', variable=self.multislice,
                                    command=self.sliceState)
        self.multislice.set(True)
        labelProj = Label(frame1, text="Number of Projections")
        self.nproj = IntVar()
        self.entryNproj = Entry(frame1, textvariable=self.nproj,
                                state=NORMAL if self.multislice.get() else DISABLED)
        self.nproj.set(900)
        labelSlice = Label(frame1, text='Slice No.')
        self.slice = IntVar()
        self.entrySlice = Entry(frame1, textvariable=self.slice,
                                state=NORMAL if self.multislice.get() else DISABLED)
        self.slice.set(0)

        # ct correction
        self.ct = BooleanVar()
        cbtCTcorrection = Checkbutton(frame1, text='Correction with Side Pixels',
                                      variable=self.ct, command=self.sideWidthState)
        self.ct.set(True)
        labelSideWidth = Label(frame1, text='Side Width (Pixels)')
        self.side_width = IntVar()
        self.entrySideWidth = Entry(frame1, textvariable=self.side_width,
                                    state=NORMAL if self.ct.get() else DISABLED)
        self.side_width.set(50)

        # energy range
        self.e_range = BooleanVar()
        cbtRange = Checkbutton(frame1, text='Energy Range', variable=self.e_range,
                               command=self.erangeState)
        self.e_range.set(False)
        self.e_low = DoubleVar()
        self.e_high = DoubleVar()
        labelElow = Label(frame1, text='Min(keV)')
        self.entryElow = Entry(frame1, textvariable=self.e_low,
                               state=NORMAL if self.e_range.get() else DISABLED)
        labelEhigh = Label(frame1, text='Max(keV)')
        self.entryEhigh = Entry(frame1, textvariable=self.e_high,
                                state=NORMAL if self.e_range.get() else DISABLED)

        # use torch
        self.use_torch = BooleanVar()
        cbtTorch = Checkbutton(frame1, text="Torch.Tensor", variable=self.use_torch)
        self.use_torch.set(True)

        # fix vertical motion
        self.fix_vertical_motion = BooleanVar()
        cbtFixVertical = Checkbutton(frame1, text='Fix Vertical Motion', variable=self.fix_vertical_motion)
        self.fix_vertical_motion.set(True)

        # snr
        self.snr = BooleanVar()
        cbtSNR = Checkbutton(frame1, text='Calculate SNR', variable=self.snr)
        self.snr.set(False)

        # algorithm
        labelAlgo = Label(frame1, text='Algorithm')
        self.algorithm = StringVar()
        rbtAlgo1 = Radiobutton(frame1, text='sKES Equation', variable=self.algorithm,
                               value='sKES_equation')
        rbtAlgo2 = Radiobutton(frame1, text='nnls', variable=self.algorithm, value='nnls')
        self.algorithm.set('sKES_equation')

        cbtUseFile.grid(row=1, column=0, sticky=W)
        cbtLowpass.grid(row=2, column=0, sticky=W)
        cbtMultislice.grid(row=3, column=0, sticky=W)
        labelProj.grid(row=3, column=1)
        self.entryNproj.grid(row=3, column=2)
        labelSlice.grid(row=3, column=3)
        self.entrySlice.grid(row=3, column=4)
        cbtCTcorrection.grid(row=4, column=0, sticky=W)
        labelSideWidth.grid(row=4, column=1, sticky=E)
        self.entrySideWidth.grid(row=4, column=2)
        cbtRange.grid(row=5, column=0, sticky=W)
        labelElow.grid(row=5, column=1, sticky=E)
        self.entryElow.grid(row=5, column=2)
        labelEhigh.grid(row=5, column=3)
        self.entryEhigh.grid(row=5, column=4)
        cbtTorch.grid(row=6, column=0, sticky=W)
        cbtFixVertical.grid(row=7, sticky=W)
        cbtSNR.grid(row=8, sticky=W)
        labelAlgo.grid(row=9, sticky=W)
        rbtAlgo1.grid(row=9, column=1, sticky=W)
        rbtAlgo2.grid(row=9, column=2, sticky=W)

        ####################### frame 1.recon  #######################
        frameRecon = Frame(frame1)
        frameRecon.grid(row=10, rowspan=3, columnspan=5, sticky=W)

        self.do_recon = BooleanVar()
        cbtRecon = Checkbutton(frameRecon, text='Reconstruction', variable=self.do_recon,
                               command=self.reconState)
        self.do_recon.set(False)
        # labelRecon = Label(frameRecon, text='Reconstruction')

        self.recon = StringVar()
        # rbtRecon0 = Radiobutton(frameRecon, text='No', variable=self.recon, value='No')
        self.rbtRecon1 = Radiobutton(frameRecon, text='SKIMAGE', variable=self.recon, value='skimage',
                                     state=NORMAL if self.do_recon.get() else DISABLED)
        self.rbtRecon2 = Radiobutton(frameRecon, text='IDL', variable=self.recon, value='idl',
                                     state=NORMAL if self.do_recon.get() else DISABLED)
        self.recon.set(None)
        self.center = IntVar()
        labelCenter = Label(frameRecon, text='Rotation Center')
        self.entryCenter = Entry(frameRecon, textvariable=self.center,
                                 state=NORMAL if self.do_recon.get() else DISABLED)
        self.center.set(0)

        # labelRecon.grid(row=0, sticky=W)
        cbtRecon.grid(row=0, column=0, sticky=W)
        self.rbtRecon1.grid(row=0, column=1, sticky=W)
        self.rbtRecon2.grid(row=0, column=2, sticky=W)
        labelCenter.grid(row=0, column=4)
        self.entryCenter.grid(row=0, column=5)

        ########################  frame 2  ###########################
        # Label(master, text='STEP 2: LOAD DATA and PRAY').grid(pady=15)
        frame2 = LabelFrame(tabKES, text = 'STEP 2: LOAD DATA and PRAY',font = ('Helvetica', '11'),
                            bd=4, relief=RIDGE,labelanchor=N)
        frame2.pack(ipadx=4)
        self.curr_path = os.getcwd()
        labelPath = Label(frame2, text="Data Directory")
        self.path = StringVar()
        entryPath = Entry(frame2, textvariable=self.path, width=75)
        btBrowse = Button(frame2, text="Browse", command=self.browseData)
        self.path.set(self.curr_path)

        labelPath.grid(sticky=W)
        entryPath.grid(row=0, column=1, padx=3)
        btBrowse.grid(row=0, column=2)

        labelSavePath = Label(frame2, text='Save Directory')
        self.save_path = StringVar()
        entrySavePath = Entry(frame2, textvariable=self.save_path, width=75)
        btBrowse2 = Button(frame2, text="Browse", command=self.browseSave)
        self.save_path.set(self.curr_path)

        labelSavePath.grid(row=1, sticky=W)
        entrySavePath.grid(row=1, column=1, padx=3)
        btBrowse2.grid(row=1, column=2)

        btRun = Button(frame2, height=2, width=30, text='MAY THE EDGE BE WITH ME',
                       bg='purple4', fg='snow', cursor='pirate', command=self.run)
        Button.config(btRun, font=('Helvetica', '11'))
        btRun.grid(row=2, columnspan=3)

        ###################  Frame3:  Reconstruction   #############################
        # Label(master, text='CT RECONSTRUCTION').grid(row=0, column=1, pady=15)
        tabRecon = ttk.Frame(master)
        frame3 = LabelFrame(tabRecon, text = 'Control',bd=4, relief=GROOVE,labelanchor=N)
        frame3.grid(row=0, column=0, sticky=NS)
        Label(frame3, text='').grid()

        frame3_1 = Frame(frame3)
        frame3_1.grid(sticky=W)
        Label(frame3_1, text='Load Sinogram', width=15).grid(row=0, sticky=W)
        btLoad = Button(frame3_1, text='Browse', command=self.browseSino).grid(row=0, column=1)

        # combobox(下拉列表) for choosing the target sinogram
        self.sino_index = '0'
        self.target_sino_id = StringVar()
        self.target_sino_id.set('0')
        self.cboxTarget = ttk.Combobox(frame3_1, width=6, textvariable=self.target_sino_id,
                                       state='disabled')
        self.cboxTarget.grid(row=0, column=2, padx=12)
        self.cboxTarget['values'] = ('0')
        self.cboxTarget.current(0)
        self.update_cbox = False  # Nothing has been loaded yet
        # Button to show selected sinogram
        self.btShowSino = Button(frame3_1, text='Show Selected', command=self.selectSino, state=DISABLED)
        self.btShowSino.grid(row=0, column=3)

        frame3_3 = Frame(frame3)
        frame3_3.grid(sticky=W)
        Label(frame3_3, text='Function', width=15).grid()
        self.reconFunc = StringVar()
        self.rbtReconSK = Radiobutton(frame3_3, text='SKIMAGE', variable=self.reconFunc, value='skimage')
        self.rbtReconIDL = Radiobutton(frame3_3, text='IDL', variable=self.reconFunc, value='idl')
        self.reconFunc.set('skimage')
        self.rbtReconSK.grid(row=0, column=1)
        self.rbtReconIDL.grid(row=0, column=2)

        frame3_4 = Frame(frame3)
        frame3_4.grid(sticky=W)
        self.rotCenter = IntVar()
        labelCenter2 = Label(frame3_4, text='Rotation Center', width=15).grid()
        self.entryCenter2 = Entry(frame3_4, textvariable=self.rotCenter)
        self.entryCenter2.grid(row=0, column=1)
        self.center.set(0)

        frame3_5 = Frame(frame3)
        frame3_5.grid(sticky=W)
        labelPixel = Label(frame3_5, text='Pixel Size(um)', width=15).grid()
        self.pixel = DoubleVar()
        entryPixel = Entry(frame3_5, textvariable=self.pixel).grid(row=0, column=1)
        rbtPixel_9 = Radiobutton(frame3_5, text='9um', variable=self.pixel, value=9).grid(row=0, column=2)
        rbtPixel_13 = Radiobutton(frame3_5, text='13.6um', variable=self.pixel, value=13.6).grid(row=0, column=3)
        self.pixel.set(13.6)

        frame3_6 = Frame(frame3)
        frame3_6.grid(sticky=W)
        Label(frame3_6, text='Save Directory', width=15).grid()
        self.recon_path = StringVar()
        frame3_6_1 = Frame(frame3)
        frame3_6_1.grid(sticky=EW)
        for x in range(10):  # configure the columns to have a non-zero weight
            Grid.columnconfigure(frame3_6_1, x, weight=1)
        entryReconPath = Entry(frame3_6_1, textvariable=self.recon_path) \
            .grid(columnspan=10, sticky=EW, padx=6)
        btReconPath = Button(frame3_6, text='Browse', command=self.browseReconPath).grid(row=0, column=1)
        self.auto_save = BooleanVar()
        self.color_auto_save = StringVar()
        self.color_auto_save.set('gray70')
        self.cbtAutoSave = Checkbutton(frame3_6, text='Auto Save', variable=self.auto_save,
                                       bg=self.color_auto_save.get(), command=self.AutoSave)
        self.cbtAutoSave.grid(row=0, column=2,padx=7)
        self.auto_save.set(False)

        frame3_7 = Frame(frame3)
        frame3_7.grid(sticky=W)
        self.display = BooleanVar()
        cbtDisplay = Checkbutton(frame3_7, text='Display Reconstructions', variable=self.display).grid()
        self.display.set(True)

        frame3_8 = Frame(frame3, width=50)
        frame3_8.grid()
        btStart = Button(frame3_8, text='START RECONSTRUCTION', command=self.runRecon)
        btStart.grid(padx=30)
        Button.config(btStart, font=('Helvetica', '11'))
        btSave = Button(frame3_8, text='SAVE', command=self.saveRecon)
        btSave.grid(row=0, column=1, padx=30, sticky=NS)
        # btSave.place()
        Button.config(btSave, font=('Helvetica', '11'))

        #####################    Text  Frames   ###############################
        print(frame1.winfo_width())
        frameT1 = Frame(tabKES, height=9, width=47)
        # frameT1.grid(row=2, column=0, sticky=S)
        frameT1.pack()
        sbText = Scrollbar(frameT1)
        self.text = Text(frameT1, height=9, width=72, yscrollcommand=sbText.set,
                         relief=SUNKEN)
        sbText.config(command=self.text.yview)
        self.text.grid(row=0, column=0, sticky=S)
        sbText.grid(row=0, column=1, sticky=NS)

        frameT2 = Frame(tabRecon, height=9, width=47)
        frameT2.grid(row=1,  column=0, sticky=S)
        sbText = Scrollbar(frameT2)
        self.text = Text(frameT2, height=9, width=47, yscrollcommand=sbText.set,
                         relief=SUNKEN)
        sbText.config(command=self.text.yview)
        self.text.grid(row=0, column=0, sticky=S)
        sbText.grid(row=0, column=1, sticky=NS)

        ####################   Canvas Frame    ##############################
        # self.loaded = BooleanVar()
        self.frameCanvas = LabelFrame(tabRecon,bd=4,relief=GROOVE,
                            text = "Display",labelanchor=N)
        self.frameCanvas.grid(row=0, rowspan=3, column=1,padx=7)

        # display an empty plot on canvas
        self.figure = Figure(figsize=(4, 4))
        figure = FigureCanvasTkAgg(self.figure, master=self.frameCanvas)
        # figure.draw()
        figure.get_tk_widget().grid()

        # ######################  Exit button   ################################
        btExit = Button(self.window, text="EXIT", fg="red", command=quit, height=1, width=4,
                        relief=FLAT, cursor='sailboat')
        btExit.grid(row=0, column=0, sticky=E+N,padx=5,pady=19)

        # photo = PhotoImage(master=canvas)
        # canvas.create_image(image=photo)

        ####################  change theme config   ##########################
        # # iterate through all the widgets
        # color = 'floral white'
        # # master.config(bg = color)
        # for widget in master.winfo_children():
        #     widget.grid(padx=5)
        #     if not widget.winfo_class() in ['Entry', 'Text']:
        #         widget.configure(bg=color, takefocus=True)
        #     for sub_widget in widget.winfo_children():
        #         if not sub_widget.winfo_class() in ['Entry', 'Text']:
        #             sub_widget.configure(bg=color, takefocus=True)
        #         for sub_sub_widget in sub_widget.winfo_children():
        #             if not sub_sub_widget.winfo_class() in ['Entry', 'Text', 'TCombobox']:
        #                 sub_sub_widget.configure(bg=color, takefocus=True)
        # btRun['bg'] = 'purple4'
        s = ttk.Style()
        s.configure('TNotebook',font='helvetica 11',padding=12)
        s.configure('TNotebook.Tab',padding=[141,5],font='helvetica 11')
        master.add(tabKES, text='Spectral Imaging')
        master.add(tabRecon, text='Reconstruction')

        master.mainloop()

    ###################  Command Functions   #################################

    def call_materials_gui(self):
        gui_get_materials()
        self.select_materials.set(True)
        # print(self.materials)

    # def select(self):
    #     filename = filedialog.askopenfilename(title='Please select file')
    #     # import os
    #     filename = os.path.basename(filename).split('.')[0]  # get the filename from the full path
    #     materials = read_materials_file(filename)
    #     write_materials_file(materials, 'last')
    #     self.materials=materials

    def call_arrangement_gui(self):
        gui_get_arrangement(root=self.window, save_path='')

    def reconState(self):
        self.rbtRecon1['state'] = NORMAL if self.do_recon.get() else DISABLED
        self.rbtRecon2['state'] = NORMAL if self.do_recon.get() else DISABLED
        self.entryCenter['state'] = NORMAL if self.do_recon.get() else DISABLED
        if self.do_recon.get():
            pass
        else:
            self.recon.set(None)

    def sliceState(self):
        self.entryNproj['state'] = NORMAL if self.multislice.get() else DISABLED
        self.entrySlice['state'] = NORMAL if self.multislice.get() else DISABLED

    def sideWidthState(self):
        self.entrySideWidth['state'] = NORMAL if self.ct.get() else DISABLED

    def erangeState(self):
        if self.e_range.get():
            self.entryElow['state'] = NORMAL
            self.entryEhigh['state'] = NORMAL
        else:
            self.entryElow['state'] = DISABLED
            self.entryEhigh['state'] = DISABLED

    def processMurho(self):
        murho = mp.murho(self.name.get(), 12.658)
        print("murho " + str(murho))

    def browseData(self):
        self.curr_path = choose_path()
        self.path.set(self.curr_path)
        self.save_path.set(Path(self.curr_path)/'Save')
        print(self.path.get())
        self.text.insert(END, self.curr_path + '\n')

    def browseSave(self):
        # self.curr_path = choose_path()
        # self.save_path.set(self.curr_path)
        path = choose_path()
        self.save_path.set(path)
        print(self.save_path.get())

    def run(self):
        def thread_run():
            if self.select_materials.get()==True:
                materials = read_materials_file('last')
            else:
                materials = read_materials_file('default')
            data_path = self.path.get()
            save_path = self.save_path.get()
            n_proj = self.nproj.get()
            algorithm = self.algorithm.get()
            multislice = self.multislice.get()
            slice = self.slice.get()
            ct = self.ct.get()
            side_width = self.side_width.get()
            if self.e_range.get():
                energy_range = [self.e_low.get(), self.e_high.get()]
            else:
                energy_range = 0
            lowpass = self.lowpass.get()
            use_torch = self.use_torch.get()
            use_file = self.use_file.get()
            fix_vertical_motion = self.fix_vertical_motion.get()
            if self.recon.get() == 'None':
                reconstruction = None
            else:
                reconstruction = self.recon.get()
            ct_center = self.center.get()
            snr = self.snr.get()
            _ = nei(materials=materials,data_path=data_path, save_path=save_path, n_proj=n_proj, algorithm=algorithm,
                    multislice=multislice, slice=slice, ct=ct, side_width=side_width,
                    e_range=energy_range, lowpass=lowpass, use_torch=use_torch, use_file=use_file,
                    fix_vertical_motion=fix_vertical_motion, reconstruction=reconstruction,
                    ct_center=ct_center, snr=snr)
            self.text.insert(END, 'NEI finished.\n')

        # multithreading. Prevent the program from frozen
        threading.Thread(target=thread_run).start()


    def browseSino(self):
        """
        Accepts data from '.pkl', numpy.ndarray, or image data
        :return:
        """
        filename = choose_file()
        if not filename == '':  # if something is choosen in the popup window for the save path of recons
            def thread_loadSino():
                self.sino_file = load_object(filename)
                self.text.insert(END, '\nLoaded:\n' + filename + '\n')
                # if self.sino_file is object("nei" structure), find the sino only.
                # Here the top module name is 'toolkit', because the 'save_object' function is in 'toolkit.py'
                if type(self.sino_file).__module__ == 'toolkit':
                    self.sino_file = self.sino_file.rho_t
                # if self.sino_file is from an image, pop up a warning for value change
                elif type(self.sino_file).__module__ == 'imageio.core.util':
                    self.sino_file = np.array(self.sino_file)
                    self.text.insert(END, '\nWarning:\nSinogram loaded from an IMAGE file. '
                                          r'Values are most likely not representing the true $\rho_t$.')
                # otherwise, only numpy.ndarray is supported
                elif not type(self.sino_file).__module__ == 'numpy':
                    self.text.insert(END, '\nError:\n'
                                          'Loaded Sinogram data type is NOT SUPPORTED.')
                    raise Exception('Loaded Sinogram data type is NOT SUPPORTED.')
                # use self.loaded as a switch for target sino and sino display
                # self.loaded.set(True)
                self.sino_file_shape = self.sino_file.shape
                self.num_sino = int(np.prod(np.array(self.sino_file_shape[:-2])))
                print('num_sino', self.num_sino)

                # use the number of dimensions before the last two to calculate the number of sinograms.
                # When the dimension is only two, np.prod(np.array([]))=1.0
                shape = self.sino_file_shape
                # The sinogram to reconstruct
                self.sino_array = self.sino_file.reshape(self.num_sino, shape[-2], shape[-1])
                # update the available values for target sinogram
                self.update_cbox = True
                self.cboxTarget.configure(state='readonly')
                self.btShowSino.configure(state=NORMAL)
                self.updateCbox()
                # display the first sinogram in loaded sinograms
                self.selectSino()
                self.text.insert(END, str(self.num_sino) + ' Sinograms Loaded\n')

            self.text.insert(END, 'Loading Sinograms\n')
            threading.Thread(target=thread_loadSino).start()
        else:  # do nothing
            pass

    def updateCbox(self):
        if self.update_cbox:  # if sinograms loaded, update the Combobox
            self.sino_index = ['All']
            self.sino_index.extend(list((np.arange(self.num_sino) + 1).astype(str)))
            print('sino_index', self.sino_index)
            self.cboxTarget.config(values=self.sino_index)
            self.cboxTarget.current(0)  # use the first element as the default in combobox
            print(self.cboxTarget.get())

    def selectSino(self):
        if self.update_cbox:  # Only when something has been loaded
            id = self.target_sino_id.get()
            if id == 'All':
                id = 1
            else:
                id = int(id)
            print('id', id)
            self.target_sino = self.sino_array[id - 1]
            #################  display self.target_sino  #################################
            # if target_id is 0, display the first one. (Recon for all)

            # prepare the figure to show
            self.figure = Figure(figsize=(4, 4))
            a = self.figure.add_subplot(1, 1, 1)
            a.imshow(self.target_sino)

            print('display target_sino')
            # clear the frameCanvas before displaying
            for child_widget in self.frameCanvas.winfo_children():
                child_widget.destroy()
            # display the self.figure on canvas
            # print('1')
            figure = FigureCanvasTkAgg(self.figure, master=self.frameCanvas)
            # print('2')
            # figure.draw() # this line was commented out for fixing an issue, which is when a sinogram should be dispalyed, the whole program would quit. Maybe caused by some updated modules.
            # print('3')
            figure.get_tk_widget().grid()
            # print('selectSino')

    def runRecon(self):

        def thread_runRecon():
            recon_funcs = {'idl': idl_recon, 'skimage': skimage_recon}
            ######### check how many recons need to be done
            targets = self.target_sino_id.get()
            if targets in ['All', '0']:
                # Run the recon. Note: convert the unit of pixel size to Centimeter
                self.recons = recon_funcs[self.reconFunc.get()](self.sino_file,
                                                                pixel_size=(self.pixel.get() / 10000),
                                                                center=self.rotCenter.get())
            else:
                # Run the recon. Note: convert the unit of pixel size to Centimeter
                self.recons = recon_funcs[self.reconFunc.get()](self.target_sino,
                                                                pixel_size=(self.pixel.get() / 10000),
                                                                center=self.rotCenter.get())
            ######### check whether to save or not
            if self.auto_save.get():
                save_recon(self.recon_path.get(), self.recons)

            ######### display reconstruction if needed
            def displayRecon(recon,counter):
                print(recon.shape)
                if recon.ndim >= 3:
                    for i in range(recon.shape[0]):
                        counter+=1  # use counter to name recon windows
                        displayRecon(recon[i],counter)
                else:

                    self.windowRecon = Toplevel()
                    self.windowRecon.title('RECON '+str(counter))
                    figure_recon = Figure(figsize=(5, 5))
                    a = figure_recon.add_subplot(1, 1, 1)
                    a.imshow(recon, cmap='gray_r')
                    # display the figure_recon on canvas
                    figure = FigureCanvasTkAgg(figure_recon, master=self.windowRecon)
                    # figure.draw()
                    figure.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
                    # display toolbar
                    toolbar = NavigationToolbar2Tk(figure, self.windowRecon)
                    toolbar.update()
                    figure._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)

                    print('end')
            if self.display.get():
                displayRecon(self.recons,0)# 0 is sent in as a counter start
            self.text.insert(END, '\nReconstruction finished\n')
        self.text.insert(END,'\nReconstruction started\n')
        threading.Thread(target=thread_runRecon).start()

    def browseReconPath(self):
        # choose path to save the reconstructed images from the tabRecon
        self.recon_path.set(choose_path())

    def AutoSave(self):
        # Change back ground color for the checkbutton
        Checkbutton.config(self.cbtAutoSave, bg='MediumPurple1' if self.auto_save.get() else 'gray70')
        # if (auto_save checked) & (path to save not defined), define the path
        if self.auto_save.get() and not self.recon_path.get():
            self.recon_path.set(choose_path())

    def saveRecon(self):
        # save when there is path and recon
        if not self.recon_path.get():  # if 'path' is empty, choose path.
            self.recon_path.set(choose_path())
        # check the path again, in case nothing was choosen from the last command.
        if self.recon_path.get() and self.recons:  # if 'recons' is not empty, then we can save recons
            save_recon(self.recon_path.get(), self.recons)

    def test(self):
        print('good')


NearEdgeImaging()
