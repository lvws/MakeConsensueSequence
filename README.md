# 序列一致性构建软件
## 功能和特色
* 常规方法构建SCS、DCS；
* 针对UMI可能的测序错误进行优化，合并判定为错误的UMI（ MisMatch 1 , 且不是该区域主要UMI ）到主要UMI簇中；
* SCS构建时优化细节，多条序列构建SCS时，去除少量含有异常DEL、INS的序列，以减少错位造成的大量错配。