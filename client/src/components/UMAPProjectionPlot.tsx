import * as React from "react";
import { ProjectionPlot, IProjectionPlotProps } from "./ProjectionPlot";
import { Button, Modal } from "react-bootstrap";
import { StructureCardGrid } from "./StructureCardGrid";
import { IParticle } from "../interfaces";
import { extent } from "d3-array";

export const UMAPProjectionPlot = ({
  umapPoints,
  colorCoding,
  ...innerProps
}: {
  umapPoints: IParticle[] | null;
  colorCoding: string | null;
} & IProjectionPlotProps) => {
  const [pointSize, setPointSize] = React.useState<number>(1.5);
  const [chemblSelection, setChemblSelection] = React.useState<IParticle[] | null>(null);
  const chemblUmapTrace = React.useMemo<Plotly.Data | null>(() => {
    if (!umapPoints) {
      return null;
    }
    const color = colorCoding ? umapPoints.map((p) => p.properties?.[colorCoding]) : undefined;
    console.log(umapPoints);
    console.log(colorCoding);
    // @ts-ignore
    let colorExtent = color ? extent(color) : undefined;
    return {
      x: (umapPoints || []).map((p) => p.projection[0]),
      y: (umapPoints || []).map((p) => p.projection[1]),
      text: (umapPoints || []).map((p) => p.structure),
      type: "scattergl",
      mode: "markers",
      name: "ChEMBL UMAP",
      hoverinfo: "none", // use skip to avoid hover
      marker: {
        color,
        cmin: colorExtent?.[0],
        cmax: colorExtent?.[1],
        size: pointSize,
        opacity: 0.6,
        colorscale: "Portland",
      },
    };
  }, [umapPoints, colorCoding, pointSize]);

  console.log(chemblUmapTrace);

  return (
    <>
      {/* <Modal
        show={Boolean(chemblSelection)}
        onHide={() => setChemblSelection(null)}
        size="xl"
        dialogClassName="modal-full-width"
      >
        {chemblSelection ? (
          <>
            <Modal.Header closeButton>
              <Modal.Title>View {chemblSelection.length} entries</Modal.Title>
            </Modal.Header>
            <Modal.Body
              style={{
                display: "flex",
                flexDirection: "column",
              }}
            >
              <StructureCardGrid structures={chemblSelection} />
            </Modal.Body>
            <Modal.Footer>
              <Button variant="secondary" onClick={() => setChemblSelection(null)}>
                Close
              </Button>
            </Modal.Footer>
          </>
        ) : null}
      </Modal>
      <form className="form-inline ml-5">
        <label className="mx-sm-3">ChEMBL Options:</label>
        <div className="form-group">
          <label htmlFor="chemblPointSizeInput">Size</label>
          <input
            type="range"
            className="form-control mx-sm-3"
            id="chemblPointSizeInput"
            min={1}
            max={10}
            value={pointSize}
            onChange={(e) => {
              setPointSize(e.currentTarget.valueAsNumber);
            }}
          />
        </div>
        <button
          className="btn btn-sm btn-light ml-5"
          title="Show ChEMBL in table view"
          onClick={(e) => {
            e.preventDefault();
            if (umapPoints) {
              const selectedPoints: number[] =
                // @ts-ignore
                chemblUmapTrace?.selectedpoints || [];
              const points = selectedPoints.map((p) => umapPoints?.[p]);
              // TODO: Make IUMAPParticle compatible
              setChemblSelection(points.length > 0 ? points : null);
            } else {
              setChemblSelection(null);
            }
          }}
        >
          <i className="fas fa-fw fa-table" />
        </button>
      </form>
      <ProjectionPlot title="ChEMBL UMAP" {...innerProps} additionalTraces={chemblUmapTrace} /> */}
      TODO
    </>
  );
};
