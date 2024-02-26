import numpy as np
cimport numpy as np

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

def redistribute_counts(np.ndarray[DTYPE_t, ndim=2] hist2d, 
                        np.ndarray[DTYPE_t, ndim=1] x_edges_old, 
                        np.ndarray[DTYPE_t, ndim=1] y_edges_old, 
                        np.ndarray[DTYPE_t, ndim=1] x_edges_new, 
                        np.ndarray[DTYPE_t, ndim=1] y_edges_new):

    cdef int i, j, new_i, new_j
    cdef int num_x_old = len(x_edges_old) - 1
    cdef int num_y_old = len(y_edges_old) - 1
    cdef int num_x_new = len(x_edges_new) - 1
    cdef int num_y_new = len(y_edges_new) - 1

    cdef np.ndarray[DTYPE_t, ndim=2] new_hist = np.zeros((num_x_new, num_y_new))

    cdef double total_area, x_overlap, y_overlap, normalized_overlap

    for i in range(num_x_old):
        for j in range(num_y_old):
            for new_i in range(num_x_new):
                for new_j in range(num_y_new):
                    x_overlap = max(0, min(x_edges_old[i+1], x_edges_new[new_i+1]) - max(x_edges_old[i], x_edges_new[new_i]))
                    y_overlap = max(0, min(y_edges_old[j+1], y_edges_new[new_j+1]) - max(y_edges_old[j], y_edges_new[new_j]))

                    total_area = (x_edges_old[i+1] - x_edges_old[i]) * (y_edges_old[j+1] - y_edges_old[j])
                    normalized_overlap = (x_overlap * y_overlap) / total_area

                    new_hist[new_i, new_j] += hist2d[i, j] * normalized_overlap

    return np.asarray(new_hist)
