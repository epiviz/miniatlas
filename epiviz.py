from epivizfileserver import setup_app, create_fileHandler, MeasurementManager
from epivizfileserver.trackhub import TrackHub
import os
import numpy
import pickle
import math

# custom function to normalize CH methylation
def norm(col):
    mean = numpy.mean(col)
    ncol = []
    for tcol in col:
        ncol.append(tcol - mean)
    print(ncol)
    return(ncol)

if __name__ == "__main__":
    # create measurements to load multiple trackhubs or configuration files
    mMgr = MeasurementManager()

    # create file handler, enables parallel processing of multiple requests
    mHandler = create_fileHandler()

    # add genome. - for supported genomes 
    # check https://obj.umiacs.umd.edu/genomes/index.html
    genome = mMgr.add_genome("mm10")

    # load ATAC trackhub
    th1 = TrackHub("http://data.nemoarchive.org/nemoHub/Ecker_ATAC_MajorCluster/")
    for m in th1.measurements:
        m.filehandler = mHandler
        mMgr.measurements.append(m)

    # selecting ATAC files by major cell types
    # TODO: this comes from a cluster file in the future
    atac_glut = ["Ecker_ATAC_MajorCluster_L2-3_IT_MajorCluster", "Ecker_ATAC_MajorCluster_L5_IT_Rspo1_MajorCluster", "Ecker_ATAC_MajorCluster_L5_IT_Rspo2_MajorCluster", "Ecker_ATAC_MajorCluster_L5_IT_S100b_MajorCluster", "Ecker_ATAC_MajorCluster_L6_IT_MajorCluster", "Ecker_ATAC_MajorCluster_L6_CT_MajorCluster", "Ecker_ATAC_MajorCluster_L6_NP_MajorCluster", "Ecker_ATAC_MajorCluster_L5_ET_MajorCluster"]
    atac_gaba = ["Ecker_ATAC_MajorCluster_Lamp5_MajorCluster", "Ecker_ATAC_MajorCluster_Sst_MajorCluster", "Ecker_ATAC_MajorCluster_Vip_MajorCluster", "Ecker_ATAC_MajorCluster_Pvalb_Calb1_MajorCluster", "Ecker_ATAC_MajorCluster_Pvalb_Reln_MajorCluster"]

    atac_glut_mea = [t for t in th1.measurements if t.mid in atac_glut]
    atac_gaba_mea = [t for t in th1.measurements if t.mid in atac_gaba]

    # load CH Methylation Tracks
    th2 = TrackHub("http://data.nemoarchive.org/nemoHub/Ecker_CG_Methylation_MajorCluster/")
    for m in th2.measurements:
        m.filehandler = mHandler
        mMgr.measurements.append(m)

    # selecting Mehtylation tracks by Major cell type
    # TODO: this comes from a cluster file in the future
    methy_glut = ["Ecker_CG_Methylation_MajorCluster_L2-3_IT_MajorCluster", "Ecker_CG_Methylation_MajorCluster_L5_IT_Rspo1_MajorCluster", "Ecker_CG_Methylation_MajorCluster_L5_IT_Rspo2_MajorCluster", "Ecker_CG_Methylation_MajorCluster_L5_IT_S100b_MajorCluster", "Ecker_CG_Methylation_MajorCluster_L6_IT_MajorCluster", "Ecker_CG_Methylation_MajorCluster_L6_CT_MajorCluster", "Ecker_CG_Methylation_MajorCluster_L6_NP_MajorCluster", "Ecker_CG_Methylation_MajorCluster_L5_ET_MajorCluster"]
    methy_gaba = ["Ecker_CG_Methylation_MajorCluster_Lamp5_MajorCluster", "Ecker_CG_Methylation_MajorCluster_Sst_MajorCluster", "Ecker_CG_Methylation_MajorCluster_Vip_MajorCluster", "Ecker_CG_Methylation_MajorCluster_Pvalb_Calb1_MajorCluster", "Ecker_CG_Methylation_MajorCluster_Pvalb_Reln_MajorCluster"]

    methy_glut_mea = [t for t in th2.measurements if t.mid in methy_glut]
    methy_gaba_mea = [t for t in th2.measurements if t.mid in methy_gaba]

    # Add computed measurements - 
    # Compute mean of ATAC and Methylation signal per cell type
    average_atac_glut = mMgr.add_computed_measurement("computed", "mean_atac_glut", "Mean ATAC Signal (Glutamatergic)", measurements=atac_glut_mea, 
                    computeFunc=numpy.mean, annotation=None)

    average_atac_gaba = mMgr.add_computed_measurement("computed", "mean_atac_gaba", "Mean ATAC Signal (GABAergic)", measurements=atac_gaba_mea, 
                    computeFunc=numpy.mean, annotation=None)

    average_methy_glut = mMgr.add_computed_measurement("computed", "mean_methy_glut", "Mean Methy Signal (Glutamatergic)", measurements=methy_glut_mea, 
                    computeFunc=numpy.mean, annotation=None)
                
    average_methy_gaba = mMgr.add_computed_measurement("computed", "mean_methy_gaba", "Mean Methy Signal (GABAergic)", measurements=methy_gaba_mea, 
                    computeFunc=numpy.mean, annotation=None)

    # load CH Methylation
    th3 = TrackHub("http://data.nemoarchive.org/nemoHub/Ecker_CH_Methylation_MajorCluster/")
    for m in th3.measurements:
        m.filehandler = mHandler
        mMgr.measurements.append(m)
        # Compute normalized CH Methylation track using the norm function
        mMgr.add_computed_measurement("computed", m.mid + "_norm", m.name + " - Norm", measurements=[m], 
                    computeFunc=norm, computeAxis=0, annotation=m.annotation)

    # Load ATAC peaks
    fmeasurements = mMgr.import_files(os.getcwd() + "/peaks.json", mHandler)

    # setup the app from the measurements manager 
    # and run the app
    app = setup_app(mMgr)
    app.run(host="0.0.0.0", port=8000)