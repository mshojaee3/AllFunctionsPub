# ============================================================
# Export_Abaqus2CSV : Robust ODB -> 3 CSV files (2D/3D, Py2/Py3)
#   1) <prefix>_Connectivity.csv
#   2) <prefix>_NodalData.csv
#   3) <prefix>_GaussData.csv
#
# Usage (inside Abaqus/CAE noGUI script):
#   exp = Export_Abaqus2CSV(JOB + '.odb', out_prefix=JOB, step_name=STEP)
#   exp.export(connectivity=True, nodal=True, gauss=True,
#              export_X=True, export_U=True, export_E=True, export_LE=True, export_S=True)
# ============================================================

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
        inst_names = list(self.odb.rootAssembly.instances.keys())
    
        if self.instance_key is None:
            prefer = []
            if self.prefer_instance_contains:
                prefer = [k for k in inst_names if self.prefer_instance_contains.upper() in k.upper()]
            self.instance_key = prefer[0] if len(prefer) else inst_names[0]
    
        self.inst = self.odb.rootAssembly.instances[self.instance_key]
    
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

    # -------------------------
    # High-level export
    # -------------------------
    def export(self, connectivity=True, nodal=True, gauss=True,
               export_X=True, export_U=True, export_E=False, export_LE=True, export_S=False,
               export_gauss_coords=True,
               nodal_tensor_position='NODAL'):
        try:
            self.open()
    
            if connectivity:
                self.export_connectivity()
    
            if nodal:
                self.export_nodal(export_X=export_X,
                                  export_U=export_U,
                                  export_E=export_E,
                                  export_LE=export_LE,
                                  export_S=export_S,
                                  nodal_tensor_position=nodal_tensor_position)
    
            if gauss:
                self.export_gauss(export_E=export_E,
                                  export_LE=export_LE,
                                  export_S=export_S,
                                  export_coords=export_gauss_coords)
        finally:
            self.close()


# -------------------------
# Example call
# -------------------------
# exp = Export_Abaqus2CSV(JOB + '.odb', out_prefix=JOB, step_name=STEP, frame_index=-1)
# exp.export(connectivity=True, nodal=True, gauss=True,
#            export_X=True, export_U=True, export_E=True, export_LE=True, export_S=True)
