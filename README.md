# PDNET
## A fully open-source framework for deep learning protein real-valued distances

As deep learning algorithms drive the progress in protein structure prediction, a lot remains to be studied at this merging superhighway of deep learning and protein structure prediction. Recent findings show that inter-residue distance prediction, a more granular version of the well-known contact prediction problem, is a key to predicting accurate models. However, deep learning methods that predict these distances are still in the early stages of their development. To advance these methods and develop other novel methods, a need exists for a small and representative dataset packaged for faster development and testing. In this work, we introduce protein distance net (PDNET), a framework that consists of one such representative dataset along with the scripts for training and testing deep learning methods. The framework also includes all the scripts that were used to curate the dataset, and generate the input features and distance maps. Deep learning models can also be trained and tested in a web browser using free platforms such as Google Colab. We discuss how the PDNET framework can be used to predict contacts, distance intervals, and real-valued distances.

## Full dataset
http://deep.cs.umsl.edu/pdnet/  

## Distance prediction compared with the image depth prediction problem
![](./depth_pred_comparison.png)
(Figure above) Comparison of the protein inter-residue distance prediction problem with the 'depth prediction from single
image problem' in computer vision. In both problems the input to the deep learning model is a volume and the
output is a 2D matrix. The depth predictions for this specific image (top right corner) were obtained by running the
pretrained FCRN method.

## Please cite
"A fully open-source framework for deep learning protein real-valued distances", B. Adhikari, Scientific Reports, 2020.
DOI: [https://www.nature.com/articles/s41598-020-70181-0#citeas](https://www.nature.com/articles/s41598-020-70181-0#citeas)

## Watch to learn more about PDNET
[https://www.youtube.com/watch?v=uAIuA1O7iE8](https://www.youtube.com/watch?v=uAIuA1O7iE8)

## Where to start?
### In a server without notebook
```bash
wget http://deep.cs.umsl.edu/pdnet/train-data.tar.gz
tar zxvf train-data.tar.gz
wget http://deep.cs.umsl.edu/pdnet/train-src.tar.gz
tar zxvf train-src.tar.gz
cd deep-learning/src/
python3 train.py -t distance -w distance.hdf5 -n 300 -c 128 -e 64 -d 16 -f 64 -p ../../data/ -v 0 -o 0 
```
### In Google Colab
Open the `pdnet_distance.ipynb` file inside the `notebooks` folder in [Google Colab](https://colab.research.google.com/) and select a GPU runtime environment. If you are new to Google Colab, please watch [this](https://www.youtube.com/watch?v=PVsS9WtwVB8).

## Contact
Badri Adhikari  
https://badriadhikari.github.io/  
