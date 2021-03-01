import * as React from "react";
import { Modal, Button } from "react-bootstrap";
import { CurveEditor, DEFAULT_CURVE_TYPES, IPoint } from "visyn_component_curve_editor";
// Include style for default styling
import "visyn_component_curve_editor/dist/curveEditor.css";
import { scaleLinear } from "d3-scale";
import { extent } from "d3-array";

const MARGIN_BOTTOM = 50;
const MARGIN_RIGHT = 50;

const VIEWBOX_SIZE: [number, number] = [600, 400];

export function CurveEditorModal({
  open,
  initialPoints = [],
  setOpen,
  onSave,
}: {
  open: boolean;
  initialPoints?: { x: number; y: number }[];
  setOpen?(show: boolean): void;
  onSave?(value: { x: number; y: number }[]): void;
}) {
  const [points, setPoints] = React.useState<IPoint[]>([]);
  const [min, setMin] = React.useState<number>(0);
  const [max, setMax] = React.useState<number>(1);
  const [inputMin, setInputMin] = React.useState<number>(min);
  const [inputMax, setInputMax] = React.useState<number>(max);

  React.useEffect(() => {
    setInputMin(min);
    setInputMax(max);
  }, [min, max]);

  React.useEffect(() => {
    setPoints(initialPoints.map(({ x, y }, i) => ({ id: i.toString(), x, y })));
    const [minimum, maximum] = extent(initialPoints.map(({ x }) => x));
    setMin(minimum ?? 0);
    setMax(maximum ?? 1);
  }, [initialPoints]);

  const scales = React.useMemo(() => {
    return {
      x: scaleLinear()
        .domain([min, max])
        .range([0, VIEWBOX_SIZE[0] - 2 * MARGIN_RIGHT])
        .clamp(true),
      y: scaleLinear()
        .domain([0, 1])
        .range([VIEWBOX_SIZE[1] - 2 * MARGIN_BOTTOM, 0])
        .clamp(true),
    };
  }, [min, max]);

  const handleClose = React.useCallback(() => {
    setOpen?.(false);
  }, [setOpen]);

  const handleSave = React.useCallback(() => {
    onSave?.(points.map(({ x, y }) => ({ x, y })));
  }, [onSave, points]);

  return (
    <>
      <Modal show={open} onHide={handleClose} size="lg">
        <Modal.Header closeButton>
          <Modal.Title>Draw desirability</Modal.Title>
        </Modal.Header>
        <Modal.Body>
          <form
            className="form-inline mb-2"
            style={{ alignItems: "baseline", flexFlow: "row", justifyContent: "center" }}
            onSubmit={(e) => {
              e.preventDefault();
              e.stopPropagation();
              setMin(inputMin);
              setMax(inputMax);
            }}
          >
            <div className="form-group mr-4">
              <label htmlFor="curveMinInput">Minimum</label>
              <input
                type="number"
                className="form-control form-control-sm ml-2"
                id="curveMinInput"
                value={inputMin}
                onChange={(e) => setInputMin(e.currentTarget.valueAsNumber)}
              />
            </div>
            <div className="form-group mr-4">
              <label htmlFor="curveMaxInput">Maximum</label>
              <input
                type="number"
                className="form-control  form-control-sm ml-2"
                id="curveMaxInput"
                value={inputMax}
                onChange={(e) => setInputMax(e.currentTarget.valueAsNumber)}
              />
            </div>
            <button type="submit" className="btn btn-sm btn-primary">
              Apply
            </button>
          </form>
          <CurveEditor
            points={points}
            setPoints={setPoints}
            scales={scales}
            viewBoxSize={VIEWBOX_SIZE}
            curveType={DEFAULT_CURVE_TYPES["Linear"]}
          />
        </Modal.Body>
        <Modal.Footer>
          <Button variant="secondary" onClick={handleClose}>
            Close
          </Button>
          <Button variant="primary" onClick={handleSave}>
            Save
          </Button>
        </Modal.Footer>
      </Modal>
    </>
  );
}
