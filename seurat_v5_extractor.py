from core.seurat_extractor import SeuratExtractor


class SeuratV5Extractor(SeuratExtractor):

    def _extract_assay(self, data, assay_name, layer_name):

        assay_pointer = data['assays']['data'][assay_name]['data']
        X = assay_pointer['layers']['data'][layer_name].T
        var_names = assay_pointer['features']['data']
        return var_names, X