# -------------------------
# Example call
# -------------------------
# RAW_BASE = "https://raw.githubusercontent.com/mshojaee3/AllFunctionsPub/main/"
# def Py_1A_download_lib_from_github(py_name):
#     # folder of this script (fallback to cwd)
#     try:
#         d = os.path.dirname(os.path.abspath(__file__))
#     except Exception:
#         d = os.getcwd()
#     # save locally with spaces removed (import-friendly)
#     local_name = py_name.replace(' ', '')
#     local_path = os.path.join(d, local_name)
#     if os.path.isfile(local_path):
#         return local_path
#     url = RAW_BASE + py_name.replace(' ', '%20')
#     try:
#         try:
#             from urllib.request import urlopen  # Py3
#         except ImportError:
#             from urllib2 import urlopen         # Py2
#         data = urlopen(url).read()
#         with open(local_path, 'wb') as f:
#             f.write(data)
#         return local_path
#     except Exception as e:
#         print("Could not download:", url, "|", str(e))
#         return None
#
#
#Py_1A_download_lib_from_github("Abaqus_1A_Export_Results.py")
#from Abaqus_1A_Export_Results import Export_Abaqus2CSV
#exp = Export_Abaqus2CSV(JOB + '.odb', out_prefix=JOB, step_name=STEP)
#exp.export(connectivity=True, nodal=True, gauss=True,
#           export_X=True, export_U=True, export_E=True, export_LE=True, export_S=True,
#           ALLSE=True,
#           frames='last')  # 'all'


from abaqusConstants import *
import csv

class Export_Abaqus2CSV(object):

    def __init__(self, odb_or_path, out_prefix='Main', step_name=None, frame_index=-1,
                 instance_key=None, prefer_instance_contains='RVE'):
        self.odb_or_path = odb_or_path
        self.out_prefix = out_prefix
        self.step_name = step_name
        self.frame_index = frame_index
        self.instance_key = instance_key
        self.prefer_instance_contains = prefer_instance_contains

        self.odb = None
        self.close_odb_on_done = False
        self.inst = None
        self.frame = None
        self.coord_dim = 3

    # -------------------------
    # Safe IO
    # -------------------------
    def open_csv_safe(self, path):
        try:
            return open(path, 'w', newline='')   # Py3
        except TypeError:
            return open(path, 'wb')             # Py2

    # -------------------------
    # Safe Abaqus access
    # -------------------------
    def _safe_get_field(self, frame, name):
        try:
            return frame.fieldOutputs[name]
        except:
            return None

    def _subset_safe(self, field, region=None, position=None):
        if field is None:
            return None
        try:
            if region is not None and position is not None:
                return field.getSubset(region=region, position=position)
            if region is not None:
                return field.getSubset(region=region)
            if position is not None:
                return field.getSubset(position=position)
            return field
        except:
            return None

    # -------------------------
    # Formatting helpers
    # -------------------------
    def _vec_labels(self, prefix, n):
        if prefix == 'X':
            base = ['X', 'Y', 'Z']
        else:
            base = [prefix+'1', prefix+'2', prefix+'3']
        return base[:n]

    def _tensor_comp_labels(self, prefix, ncomp):
        if ncomp == 3:
            return [prefix+'11', prefix+'22', prefix+'12']
        if ncomp == 4:
            return [prefix+'11', prefix+'22', prefix+'33', prefix+'12']
        if ncomp == 6:
            return [prefix+'11', prefix+'22', prefix+'33', prefix+'12', prefix+'13', prefix+'23']
        return [prefix+str(i+1) for i in range(ncomp)]

    # -------------------------
    # ODB open/close + selection
    # -------------------------
    def open(self):
        # Accept opened odb object or a path
        if hasattr(self.odb_or_path, 'steps') and hasattr(self.odb_or_path, 'rootAssembly'):
            self.odb = self.odb_or_path
            self.close_odb_on_done = False
        else:
            # Try session.openOdb first (CAE), then odbAccess.openOdb (noGUI safe)
            odb_obj = None
            try:
                from abaqus import session  # works in CAE/noGUI if abaqus module provides session
                odb_obj = session.openOdb(self.odb_or_path)
            except:
                try:
                    from odbAccess import openOdb
                    odb_obj = openOdb(path=self.odb_or_path, readOnly=False)
                except Exception as e:
                    raise RuntimeError("Could not open ODB: %s  (%s)" % (self.odb_or_path, str(e)))
    
            self.odb = odb_obj
            self.close_odb_on_done = True
    
        # Instance selection
        # Instance selection (ROBUST)
        inst_names = list(self.odb.rootAssembly.instances.keys())
        inst_objs  = self.odb.rootAssembly.instances
        
        # If user explicitly set instance_key, use it
        if self.instance_key is not None:
            if self.instance_key not in inst_objs:
                raise RuntimeError("instance_key '%s' not found. Available: %s"
                                   % (self.instance_key, inst_names))
            self.inst = inst_objs[self.instance_key]
        else:
            # Optional: filter by substring first (like your old BRICK logic)
            candidates = inst_names
            if self.prefer_instance_contains:
                s = self.prefer_instance_contains.upper()
                filtered = [k for k in inst_names if s in k.upper()]
                if filtered:
                    candidates = filtered
        
            # Pick instance with max elements, tie-break by nodes
            best = None
            best_key = None
            for k in candidates:
                inst = inst_objs[k]
                try:
                    ne = len(inst.elements)
                except:
                    ne = 0
                try:
                    nn = len(inst.nodes)
                except:
                    nn = 0
                score = (ne, nn)
                if (best is None) or (score > best):
                    best = score
                    best_key = k
        
            if best_key is None:
                raise RuntimeError("No instances found in ODB.")
            self.instance_key = best_key
            self.inst = inst_objs[best_key]
        
        print("Using instance:", self.inst.name, "nodes:", len(self.inst.nodes), "elements:", len(self.inst.elements))
    
    
        # Step + frame selection
        step_keys = list(self.odb.steps.keys())
        if self.step_name is None:
            self.step_name = step_keys[-1]
        step = self.odb.steps[self.step_name]
        self.frame = step.frames[self.frame_index]
    
        # Coord dimension
        try:
            self.coord_dim = len(self.inst.nodes[0].coordinates)
        except:
            self.coord_dim = 3
    
        print('--- Export_Abaqus2CSV OPEN ---')
        print(' ODB:', self.odb.name)
        print(' Instance:', self.instance_key)
        print(' Step:', self.step_name, 'Frame:', self.frame_index)
        print(' CoordDim:', self.coord_dim)

    def close(self):
        if self.odb is not None and self.close_odb_on_done:
            try:
                self.odb.close()
            except:
                pass
        self.odb = None
        self.inst = None
        self.frame = None

    # -------------------------
    # Export: Connectivity
    # -------------------------
    def export_connectivity(self, filename=None):
        if filename is None:
            filename = self.out_prefix + '_Connectivity.csv'
    
        # numeric TypeID mapping (same as your old script)
        type_map = {
            'CPS3': 3, 'CPE3': 3,
            'CPS6': 6, 'CPE6': 6,
            'CPS4': 4, 'CPE4': 4,
            'CPS8': 8, 'CPE8': 8,
            # add more if you use them:
            'C3D4': 4, 'C3D10': 10, 'C3D8': 8, 'C3D20': 20
        }
    
        # find max connectivity length for padding
        max_conn = 0
        for el in self.inst.elements:
            try:
                max_conn = max(max_conn, len(el.connectivity))
            except:
                pass
        if max_conn < 1:
            max_conn = 8
    
        with self.open_csv_safe(filename) as f:
            w = csv.writer(f)
            header = ['Element', 'TypeID'] + ['N%d' % (i+1) for i in range(max_conn)]
            w.writerow(header)
    
            for el in self.inst.elements:
                t_id = type_map.get(el.type, 0)
                conn = list(el.connectivity)
                padded = conn + [float('nan')] * (max_conn - len(conn))
                w.writerow([el.label, t_id] + padded)
    
        print(' -> Connectivity:', filename)


    # -------------------------
    # Export: Nodal Data (X + U)
    # -------------------------
    def export_nodal(self, filename=None, export_X=True, export_U=True,
                     export_E=False, export_LE=False, export_S=False,
                     nodal_tensor_position='NODAL'):
        """
        Writes one row per node.
        Exports:
          - X coords
          - U at NODAL
          - optional E/LE/S at NODAL if available,
            else uses ELEMENT_NODAL and averages across contributing elements.
        """
    
        if filename is None:
            filename = self.out_prefix + '_NodalData.csv'
    
        # -------------------------
        # Displacement (U) at nodes
        # -------------------------
        u_dict = {}
        if export_U:
            fo_u = self._subset_safe(self._safe_get_field(self.frame, 'U'),
                                     region=self.inst, position=NODAL)
            if fo_u is None:
                fo_u = self._subset_safe(self._safe_get_field(self.frame, 'U'), region=self.inst)
            if fo_u is not None:
                for v in fo_u.values:
                    u_dict[v.nodeLabel] = v.data
    
        # -------------------------
        # Helper to get nodal tensor (E/LE/S)
        # -------------------------
        def get_nodal_tensor(fname):
            """
            Returns (tensor_dict, ncomp) where tensor_dict[nodeLabel] = tuple comps
            If NODAL not available, averages ELEMENT_NODAL.
            """
            field = self._safe_get_field(self.frame, fname)
            if field is None:
                return {}, 0
    
            # Try requested position first (usually NODAL)
            pos = NODAL if str(nodal_tensor_position).upper() == 'NODAL' else ELEMENT_NODAL
            fo = self._subset_safe(field, region=self.inst, position=pos)
    
            # If not available, fallback to ELEMENT_NODAL
            used_element_nodal = False
            if fo is None or len(fo.values) == 0:
                fo = self._subset_safe(field, region=self.inst, position=ELEMENT_NODAL)
                used_element_nodal = True
    
            if fo is None or len(fo.values) == 0:
                return {}, 0
    
            # Determine number of components
            ncomp = len(fo.values[0].data)
    
            if not used_element_nodal:
                # one value per node (usually)
                d = {}
                for v in fo.values:
                    d[v.nodeLabel] = v.data
                return d, ncomp
    
            # ELEMENT_NODAL: multiple values per node -> average
            acc = {}
            cnt = {}
            for v in fo.values:
                nl = v.nodeLabel
                dat = v.data
                if nl not in acc:
                    acc[nl] = [0.0] * ncomp
                    cnt[nl] = 0
                for i in range(ncomp):
                    acc[nl][i] += dat[i]
                cnt[nl] += 1
    
            d = {}
            for nl in acc:
                c = float(cnt[nl]) if cnt[nl] else 1.0
                d[nl] = tuple([acc[nl][i] / c for i in range(ncomp)])
            return d, ncomp
    
        # -------------------------
        # Collect requested tensors
        # -------------------------
        tensors = []  # list of (fieldName, prefix, dict, ncomp)
        if export_E:
            d, n = get_nodal_tensor('E')
            if n > 0: tensors.append(('E', 'E', d, n))
        if export_LE:
            d, n = get_nodal_tensor('LE')
            if n > 0: tensors.append(('LE', 'LE', d, n))
        if export_S:
            d, n = get_nodal_tensor('S')
            if n > 0: tensors.append(('S', 'S', d, n))
    
        # -------------------------
        # Write CSV
        # -------------------------
        with self.open_csv_safe(filename) as f:
            w = csv.writer(f)
    
            header = ['Node_Label']
            if export_X:
                header += self._vec_labels('X', self.coord_dim)
            if export_U:
                header += (['U_U1','U_U2','U_U3'] if self.coord_dim == 3 else ['U_U1','U_U2'])
    
            # tensor headers
            for (_, pref, _, ncomp) in tensors:
                header += self._tensor_comp_labels(pref, ncomp)
    
            w.writerow(header)
    
            for nd in self.inst.nodes:
                lab = nd.label
                row = [lab]
    
                if export_X:
                    c = nd.coordinates
                    row += [c[i] if i < len(c) else float('nan') for i in range(self.coord_dim)]
    
                if export_U:
                    u = u_dict.get(lab, None)
                    if u is None:
                        row += ([float('nan'), float('nan'), float('nan')] if self.coord_dim == 3
                                else [float('nan'), float('nan')])
                    else:
                        if self.coord_dim == 3:
                            row += [u[0], u[1], u[2] if len(u) > 2 else 0.0]
                        else:
                            row += [u[0], u[1]]
    
                # tensors
                for (_, _, dct, ncomp) in tensors:
                    val = dct.get(lab, None)
                    if val is None:
                        row += [float('nan')] * ncomp
                    else:
                        row += [val[i] for i in range(ncomp)]
    
                w.writerow(row)
    
        print(' -> NodalData:', filename)



    # -------------------------
    # Export: Gauss Data (E, LE, S + optional coords)
    # -------------------------
    def export_gauss(self, filename=None, export_E=False, export_LE=True, export_S=False,
                     export_coords=True):
        if filename is None:
            filename = self.out_prefix + '_GaussData.csv'

        req = []
        if export_E:  req.append(('E',  'E'))
        if export_LE: req.append(('LE', 'LE'))
        if export_S:  req.append(('S',  'S'))

        # choose base field for iteration (first available tensor)
        base_values = None
        for fname, _ in req:
            fo = self._subset_safe(self._safe_get_field(self.frame, fname),
                                   region=self.inst, position=INTEGRATION_POINT)
            if fo is not None and len(fo.values) > 0:
                base_values = fo.values
                break

        # try coords at IP if requested
        coord_lookup = {}
        if export_coords:
            fo_c = self._subset_safe(self._safe_get_field(self.frame, 'COORD'),
                                     region=self.inst, position=INTEGRATION_POINT)
            if fo_c is not None and len(fo_c.values) > 0:
                for v in fo_c.values:
                    key = (v.elementLabel, v.integrationPoint)
                    cd = v.data
                    if len(cd) == 2:
                        coord_lookup[key] = (cd[0], cd[1], float('nan'))
                    else:
                        coord_lookup[key] = (cd[0], cd[1], cd[2] if len(cd) > 2 else float('nan'))

        # fallback base = coords if no tensor exists
        if base_values is None and export_coords and len(coord_lookup) > 0:
            # fabricate base from coord keys (stable order)
            base_values = []
            for (el, ip) in sorted(coord_lookup.keys()):
                class _Tmp(object): pass
                t = _Tmp()
                t.elementLabel = el
                t.integrationPoint = ip
                base_values.append(t)

        # if still no base, write header-only and return
        if base_values is None:
            with self.open_csv_safe(filename) as f:
                w = csv.writer(f)
                hdr = ['Element', 'IP']
                if export_coords:
                    hdr += ['X','Y','Z']
                w.writerow(hdr)
            print(' -> GaussData (empty):', filename)
            return

        # build (el,ip)->data dict for each requested field
        field_ip = {}
        comp_counts = {}
        headers = []

        for fname, pref in req:
            fo = self._subset_safe(self._safe_get_field(self.frame, fname),
                                   region=self.inst, position=INTEGRATION_POINT)
            if fo is None or len(fo.values) == 0:
                continue

            d = {}
            for v in fo.values:
                d[(v.elementLabel, v.integrationPoint)] = v.data
            field_ip[fname] = d

            # determine ncomp from first value
            sample = None
            for k in d:
                sample = d[k]
                break
            if sample is None:
                continue
            ncomp = len(sample)
            comp_counts[fname] = ncomp
            headers += self._tensor_comp_labels(pref, ncomp)

        with self.open_csv_safe(filename) as f:
            w = csv.writer(f)
            hdr = ['Element', 'IP']
            if export_coords:
                hdr += ['X','Y','Z']
            hdr += headers
            w.writerow(hdr)

            for v in base_values:
                el = v.elementLabel
                ip = v.integrationPoint
                key = (el, ip)

                row = [el, ip]

                if export_coords:
                    xyz = coord_lookup.get(key, (float('nan'), float('nan'), float('nan')))
                    row += [xyz[0], xyz[1], xyz[2]]

                for fname, _ in req:
                    if fname not in field_ip:
                        continue
                    ncomp = comp_counts.get(fname, 0)
                    val = field_ip[fname].get(key, None)
                    if val is None:
                        row += [float('nan')] * ncomp
                    else:
                        row += [val[i] for i in range(len(val))]

                w.writerow(row)

        print(' -> GaussData:', filename)
    

    def export_ALLSE(self, filename=None):
        """
        Export ALLSE (strain energy) history output to CSV.
        Writes: Time, ALLSE, HistoryRegion
        """
        if filename is None:
            filename = self.out_prefix + '_ALLSE.csv'
        step = self.odb.steps[self.step_name]
        rows = []
        for hr_name, hr in step.historyRegions.items():
            try:
                ho = hr.historyOutputs
            except:
                continue
            if 'ALLSE' not in ho:
                continue
            data = ho['ALLSE'].data  # [(time, value), ...]
            for (t, val) in data:
                rows.append((t, val, hr_name))
        with self.open_csv_safe(filename) as f:
            w = csv.writer(f)
            w.writerow(['Time', 'ALLSE', 'HistoryRegion'])
            for r in rows:
                w.writerow(r)
        if len(rows) == 0:
            print(" -> ALLSE not found for step '%s' (header written only)" % self.step_name)
        else:
            print(" -> ALLSE:", filename, "(%d rows)" % len(rows))
    
    # -------------------------
    # High-level export
    # -------------------------
    
    def _get_step_and_frame_indices(self, frames='last'):
        """
        frames:
          - 'last'  -> export only last frame (or self.frame_index if explicitly set)
          - 'all'   -> export all frames in step
          - list/tuple of ints -> export those frame indices
        returns: (step_obj, frame_indices)
        """
        step_keys = list(self.odb.steps.keys())
        if self.step_name is None:
            self.step_name = step_keys[-1]

        step = self.odb.steps[self.step_name]

        # Decide frame indices
        if frames == 'all':
            idxs = list(range(len(step.frames)))
        elif frames == 'last' or frames is None:
            # respect explicit self.frame_index if user set it, otherwise last
            fi = self.frame_index if self.frame_index is not None else -1
            # normalize negative
            if fi < 0:
                fi = len(step.frames) + fi
            idxs = [fi]
        else:
            # user passed explicit indices
            idxs = list(frames)
        return step, idxs
    
    
    # -------------------------
    # High-level export
    # -------------------------
    def export(self, connectivity=True, nodal=True, gauss=True,
               export_X=True, export_U=True, export_E=False, export_LE=True, export_S=False,
               export_gauss_coords=True,
               nodal_tensor_position='NODAL',
               frames='last',                 # <-- NEW: 'last' or 'all' or [0,5,9]
               suffix_with_step_frame=True,
               ALLSE=True):  
        try:
            self.open()

            # Connectivity once
            if connectivity:
                self.export_connectivity()
            
            if ALLSE:
                self.export_ALLSE(filename=self.out_prefix + '_ALLSE.csv')
            
            # Decide which frames to export
            step, frame_indices = self._get_step_and_frame_indices(frames=frames)

            for fi in frame_indices:
                self.frame_index = fi
                self.frame = step.frames[fi]

                # Filename suffix
                suf = ''
                if suffix_with_step_frame and (len(frame_indices) > 1):
                    suf = '_%s_f%04d' % (self.step_name, fi)

                if nodal:
                    self.export_nodal(
                        filename=self.out_prefix + '_NodalData' + suf + '.csv',
                        export_X=export_X,
                        export_U=export_U,
                        export_E=export_E,
                        export_LE=export_LE,
                        export_S=export_S,
                        nodal_tensor_position=nodal_tensor_position
                    )

                if gauss:
                    self.export_gauss(
                        filename=self.out_prefix + '_GaussData' + suf + '.csv',
                        export_E=export_E,
                        export_LE=export_LE,
                        export_S=export_S,
                        export_coords=export_gauss_coords
                    )

            print("Done. Exported %d frame(s) from step '%s'." % (len(frame_indices), self.step_name))

        finally:
            self.close()
