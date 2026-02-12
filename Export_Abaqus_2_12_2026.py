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
            self.odb = session.openOdb(self.odb_or_path)
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
            header = ['Element', 'ElemType'] + ['N%d' % (i+1) for i in range(max_conn)]
            w.writerow(header)

            for el in self.inst.elements:
                conn = list(el.connectivity)
                padded = conn + [float('nan')] * (max_conn - len(conn))
                w.writerow([el.label, el.type] + padded)

        print(' -> Connectivity:', filename)

    # -------------------------
    # Export: Nodal Data (X + U)
    # -------------------------
    def export_nodal(self, filename=None, export_X=True, export_U=True):
        if filename is None:
            filename = self.out_prefix + '_NodalData.csv'

        u_dict = {}
        if export_U:
            fo_u = self._subset_safe(self._safe_get_field(self.frame, 'U'),
                                     region=self.inst, position=NODAL)
            if fo_u is None:
                fo_u = self._subset_safe(self._safe_get_field(self.frame, 'U'), region=self.inst)
            if fo_u is not None:
                for v in fo_u.values:
                    u_dict[v.nodeLabel] = v.data

        with self.open_csv_safe(filename) as f:
            w = csv.writer(f)
            header = ['Node_Label']
            if export_X:
                header += self._vec_labels('X', self.coord_dim)
            if export_U:
                header += (['U1','U2','U3'] if self.coord_dim == 3 else ['U1','U2'])
            w.writerow(header)

            for nd in self.inst.nodes:
                row = [nd.label]
                c = nd.coordinates
                if export_X:
                    row += [c[i] if i < len(c) else float('nan') for i in range(self.coord_dim)]
                if export_U:
                    u = u_dict.get(nd.label, None)
                    if u is None:
                        row += ([float('nan'), float('nan'), float('nan')] if self.coord_dim == 3
                                else [float('nan'), float('nan')])
                    else:
                        if self.coord_dim == 3:
                            row += [u[0], u[1], u[2] if len(u) > 2 else 0.0]
                        else:
                            row += [u[0], u[1]]
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
               export_gauss_coords=True):
        try:
            self.open()
            if connectivity:
                self.export_connectivity()
            if nodal:
                self.export_nodal(export_X=export_X, export_U=export_U)
            if gauss:
                self.export_gauss(export_E=export_E, export_LE=export_LE,
                                  export_S=export_S, export_coords=export_gauss_coords)
        finally:
            self.close()


# -------------------------
# Example call
# -------------------------
# exp = Export_Abaqus2CSV(JOB + '.odb', out_prefix=JOB, step_name=STEP, frame_index=-1)
# exp.export(connectivity=True, nodal=True, gauss=True,
#            export_X=True, export_U=True, export_E=True, export_LE=True, export_S=True)
